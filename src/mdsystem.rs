use variables::Variables;
use parameter::Parameter;
use observer::Observer;
use utils::f64vec3::F64vec3;
use configmaker::ConfigMaker;

pub struct MDSystem {
    vars: Variables,
    param: Parameter,
    obs: Observer,
}

impl MDSystem {
    pub fn new() -> MDSystem {
        MDSystem {
            vars: Variables::new(),
            param: Parameter::new(),
            obs: Observer::new(),
        }
    }

    fn calculate_force(&mut self) {
        let vars  = &mut self.vars;
        let param = &self.param;

        let n_atoms = vars.number_of_atoms();
        let cl2 = param.cl * param.cl;
        let dt  = param.dt;
        for i in 0..n_atoms-1 {
            for j in i+1..n_atoms {
                let mut dr = F64vec3::new(
                    vars.atoms[j].r.x - vars.atoms[i].r.x,
                    vars.atoms[j].r.y - vars.atoms[i].r.y,
                    vars.atoms[j].r.z - vars.atoms[i].r.z,
                );
                param.adjust_pbc(&mut dr);
                let r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
                if r2 > cl2 {continue;}
                let r6 = r2 * r2 * r2;
                let df = (24.0 * r6 - 48.0) / (r6 * r6 * r2) * dt;
                vars.atoms[i].v.x += df * dr.x;
                vars.atoms[i].v.y += df * dr.y;
                vars.atoms[i].v.z += df * dr.z;
                vars.atoms[j].v.x -= df * dr.x;
                vars.atoms[j].v.y -= df * dr.y;
                vars.atoms[j].v.z -= df * dr.z;
            }
        }
    }

    fn update_position(&mut self) {
        let vars  = &mut self.vars;
        let param = &self.param;
        let dt_h = param.dt * 0.5;
        for atom in &mut vars.atoms {
            atom.r.x += atom.v.x * dt_h;
            atom.r.y += atom.v.y * dt_h;
            atom.r.z += atom.v.z * dt_h;
        }
    }

    fn periodic(&mut self) {
        let vars  = &mut self.vars;
        let param = &self.param;
        let bl = &param.box_len;
        for a in &mut vars.atoms {
            if a.r.x < 0.0  {a.r.x += bl.x;};
            if a.r.y < 0.0  {a.r.y += bl.y;};
            if a.r.z < 0.0  {a.r.z += bl.z;};
            if a.r.x > bl.x {a.r.x -= bl.x;};
            if a.r.y > bl.y {a.r.y -= bl.y;};
            if a.r.z > bl.z {a.r.z -= bl.z;};
        }
    }

    fn calculate(&mut self) {
        self.update_position();
        self.calculate_force();
        self.update_position();
        self.periodic();
        self.vars.time += self.param.dt;
    }

    pub fn run(&mut self) {
        ConfigMaker::make_fcc(&mut self.vars,
                              &self.param);
        self.vars.set_initial_velocity(1.0);
        let steps = 10000;
        let observe = 100;
        for i in 0..steps {
            if (i % observe) == 0 {
                self.obs.traject(&self.vars);
                let k = self.obs.kinetic_energy(&self.vars);
                let u = self.obs.potential_energy(&self.vars,
                                                  &self.param);
                println!("t = {}, k = {}, u = {}, tot = {}",
                         self.vars.time, k, u, k + u);

            }
            self.calculate()
        }
    }
}



