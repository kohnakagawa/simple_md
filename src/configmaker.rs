extern crate rand;

use variables::Variables;
use parameter::Parameter;
use utils::f64vec3::F64vec3;
use self::rand::distributions::{IndependentSample, Range};

pub trait ConfigMaker {
    fn make_conf(&self, vars: &mut Variables, param: &Parameter);
}

#[allow(dead_code)]
pub struct RandomConfigMaker {
    density: f64,
}

impl RandomConfigMaker {
    #[allow(dead_code)]
    pub fn new(density: f64) -> RandomConfigMaker {
        RandomConfigMaker {
            density: density
        }
    }
}

impl ConfigMaker for RandomConfigMaker {
    fn make_conf(&self, vars: &mut Variables, param: &Parameter) {
        let mut rng = rand::thread_rng();
        let bl = &param.box_len;
        let between_x = Range::new(0.0, bl.x);
        let between_y = Range::new(0.0, bl.y);
        let between_z = Range::new(0.0, bl.z);

        let num_atoms = (self.density * bl.x * bl.y * bl.z) as usize;
        for _ in 0..num_atoms {
            vars.add_atoms(F64vec3::new(
                between_x.ind_sample(&mut rng),
                between_y.ind_sample(&mut rng),
                between_z.ind_sample(&mut rng),
            ));
        }
    }
}

#[allow(dead_code)]
pub struct FccConfigMaker {
    density: f64,
}

impl FccConfigMaker {
    #[allow(dead_code)]
    pub fn new(density: f64) -> FccConfigMaker {
        FccConfigMaker {
            density: density
        }
    }
}

impl ConfigMaker for FccConfigMaker {
    fn make_conf(&self, vars: &mut Variables, param: &Parameter) {
        let s  = 1.0 / (self.density * 0.25).powf(1.0 / 3.0);
        let hs = s * 0.5;
        let bl = &param.box_len;
        let l = if bl.x < bl.y {bl.x} else {bl.y};
        let l = if l < bl.z {l} else {bl.z};
        let isx = (l / s) as u32;
        let isy = isx;
        let isz = isx;
        for iz in 0..isz {
            for iy in 0..isy {
                for ix in 0..isx {
                    let x = (ix as f64) * s;
                    let y = (iy as f64) * s;
                    let z = (iz as f64) * s;
                    vars.add_atoms(F64vec3::new(x,      y,      z     ));
                    vars.add_atoms(F64vec3::new(x,      y + hs, z + hs));
                    vars.add_atoms(F64vec3::new(x + hs, y,      z + hs));
                    vars.add_atoms(F64vec3::new(x + hs, y + hs, z     ));
                }
            }
        }
    }
}
