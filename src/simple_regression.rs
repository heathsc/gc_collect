/// Simple (one predictor) linear regression
use stat_functions::students_t::StudentsT;

#[derive(Debug, Copy, Clone)]
pub struct Coefficient {
    estimate: f64,
    standard_error: f64,
    df: usize,
}

impl Coefficient {
    pub fn estimate(&self) -> f64 {
        self.estimate
    }

    pub fn t_statistic(&self) -> f64 {
        self.estimate / self.standard_error
    }
    pub fn p(&self) -> f64 {
        let v = self.df as f64;
        let s = StudentsT::new(v).expect("Invalid df");
        let t = self.t_statistic();
        // 2-sided t-test
        2.0 * s.pt(-t.abs())
    }
}
#[derive(Debug)]
#[allow(unused)]
pub struct SimpleRegression {
    intercept: Coefficient,
    slope: Coefficient,
    residual_ss: f64,
    residual_df: usize,
}

impl SimpleRegression {
    pub fn slope(&self) -> &Coefficient {
        &self.slope
    }
}

#[derive(Default, Copy, Clone)]
struct RegSums {
    sum_x: f64,
    sum_x2: f64,
    sum_y: f64,
    sum_xy: f64,
}

impl RegSums {
    fn add_obs(&mut self, x: f64, y: f64) {
        self.sum_x += x;
        self.sum_x2 += x * x;
        self.sum_y += y;
        self.sum_xy += x * y;
    }
}

fn get_reg_sums(obs: &[(f64, f64)]) -> RegSums {
    let mut rsums = RegSums::default();
    for (x, y) in obs {
        rsums.add_obs(*x, *y)
    }
    rsums
}

pub fn simple_regression(obs: &[(f64, f64)]) -> anyhow::Result<SimpleRegression> {
    if obs.len() < 3 {
        Err(anyhow!(
            "Cannot obtain meaningful regression estimates with <3 observations"
        ))
    } else {
        // Assumulate sums
        let rs = get_reg_sums(obs);

        // Calculate determinant of X'X
        let n = obs.len() as f64;
        let det = n * rs.sum_x2 - rs.sum_x.powi(2);
        if det <= 0.0 {
            return Err(anyhow!("Numerical error during regression calculations"));
        }

        // Calculate regression coefficients
        let b0 = (rs.sum_x2 * rs.sum_y - rs.sum_x * rs.sum_xy) / det;
        let b1 = (n * rs.sum_xy - rs.sum_x * rs.sum_y) / det;

        // Calculate residual sum of squares
        let residual_ss = obs
            .iter()
            .map(|(x, y)| (y - b0 - x * b1).powi(2))
            .sum::<f64>();

        // Residual variance
        let res_var = residual_ss / (n - 2.0);
        let df = obs.len() - 2;

        let intercept = Coefficient {
            estimate: b0,
            standard_error: (rs.sum_x2 * res_var / det).sqrt(),
            df,
        };

        let slope = Coefficient {
            estimate: b1,
            standard_error: (n * res_var / det).sqrt(),
            df,
        };

        Ok(SimpleRegression {
            intercept,
            slope,
            residual_ss,
            residual_df: df,
        })
    }
}

mod test {
    #[allow(unused_imports)]
    use super::*;

    #[test]
    fn regression_test() {
        let obs = [(1.0, 4.0), (2.0, 6.0), (4.0, 7.0), (3.0, 6.5), (5.0, 10.0)];
        let reg = simple_regression(&obs).expect("Error in regression");
        println!("{:?}", reg);
        println!("{}", reg.slope().p());
        assert!((reg.slope().p() - 0.0140732510).abs() < 1.0e-8);
    }
}
