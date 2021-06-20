use std::f64::consts::PI;

use once_cell::unsync::Lazy;

pub const SPH_RESTDENSITY: f64 = 600.0;
pub const SPH_INTSTIFF: f64 = 1.0;
pub const SPH_PMASS: f64 = 0.00020543;
pub const SPH_SIMSCALE: f64 = 0.004;
pub const H: f64 = 0.01;
pub const DT: f64 = 0.004;
pub const SPH_VISC: f64 = 0.2;
pub const SPH_LIMIT: f64 = 200.0;
pub const SPH_RADIUS: f64 = 0.004;
pub const SPH_EPSILON: f64 = 0.00001;
pub const SPH_EXTSTIFF: f64 = 10000.0;
pub const SPH_EXTDAMP: f64 = 256.0;
pub const SPH_PDIST: Lazy<f64> = Lazy::new(|| (SPH_PMASS / SPH_RESTDENSITY).powf(1.0 / 3.0));
pub const MIN: [f64; 2] = [0.0, 0.0];
pub const MAX: [f64; 2] = [20.0, 50.0];
pub const INITMIN: [f64; 2] = [0.0, 0.0];
pub const INITMAX: [f64; 2] = [10.0, 20.0];
pub const POLY6KERN: Lazy<f64> = Lazy::new(|| 315.0 / (64.0 * PI * H.powf(9.0)));
pub const SPIKYKERN: Lazy<f64> = Lazy::new(|| -45.0 / (PI * H.powf(6.0)));
pub const LAPKERN: Lazy<f64> = Lazy::new(|| 45.0 / (PI * H.powf(6.0)));
pub const D: Lazy<f64> = Lazy::new(|| *SPH_PDIST / SPH_SIMSCALE * 0.95);
