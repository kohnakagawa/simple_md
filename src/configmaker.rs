use variables::Variables;
use parameter::Parameter;
use utils::f64vec3::F64vec3;

pub struct ConfigMaker {}

impl ConfigMaker {
    pub fn make_fcc(vars: &mut Variables, param: &Parameter) {
        let density = 0.5_f64;
        let s  = 1.0 / (density * 0.25).powf(1.0 / 3.0);
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
