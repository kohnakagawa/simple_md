use utils::f64vec3::F64vec3;

pub struct Parameter {
    pub box_len: F64vec3,
    pub dt: f64,
    pub cl: f64,
    pub c0: f64,
    pub margin: f64,
}

impl Parameter {
    pub fn new() -> Parameter {
        let cl   = 2.0;
        let rc2  = 1.0 / (cl * cl);
        let rc6  = rc2 * rc2 * rc2;
        let rc12 = rc6 * rc6;
        let c0   = -4.0 * (rc12 - rc6);
        Parameter {
            box_len: F64vec3::new(10.0, 10.0, 10.0),
            dt: 0.01,
            cl: cl,
            c0: c0,
            margin: 0.5,
        }
    }

    pub fn adjust_pbc(&self, dr: &mut F64vec3) {
        let bl = &self.box_len;
        let bh = F64vec3 {
            x: bl.x * 0.5,
            y: bl.y * 0.5,
            z: bl.z * 0.5,
        };
        if dr.x < -bh.x {dr.x += bl.x;}
        if dr.x >  bh.x {dr.x -= bl.x;}
        if dr.y < -bh.y {dr.y += bl.y;}
        if dr.y >  bh.y {dr.y -= bl.y;}
        if dr.z < -bh.z {dr.z += bl.z;}
        if dr.z >  bh.z {dr.z -= bl.z;}
    }
}
