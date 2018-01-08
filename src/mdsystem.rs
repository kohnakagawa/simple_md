use variables::Variables;
use parameter::Parameter;
use fcalculator::ForceCalculator;
use meshlist::MeshList;
use observer::Observer;
use configmaker::ConfigMaker;
use configmaker::FccConfigMaker;

pub struct MDSystem {
    vars: Variables,
    param: Parameter,
    obs: Observer,
    mesh: MeshList,
    margin_length: f64,
}

impl MDSystem {
    pub fn new() -> MDSystem {
        MDSystem {
            vars: Variables::new(),
            param: Parameter::new(),
            obs: Observer::new(),
            mesh: MeshList::new(),
            margin_length: 0.0,
        }
    }

    fn check_pairlist(&mut self) {
        let vars = &mut self.vars;
        let param = &self.param;

        let mut vmax2 = 0.0;
        for atom in &vars.atoms {
            let v2 = atom.v.norm2();
            if vmax2 < v2 {vmax2 = v2}
        }
        let vmax = vmax2.sqrt();

        self.margin_length -= vmax * param.dt * 2.0;
        if self.margin_length < 0.0 {
            self.margin_length = param.margin;
            self.mesh.make_sorted_list(vars, &self.param);
        }
    }

    fn calculate_force(&mut self) {
        ForceCalculator::calculate_force(&mut self.vars,
                                         &self.mesh,
                                         &self.param);
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
        self.check_pairlist();
        self.calculate_force();
        self.update_position();
        self.periodic();
        self.vars.time += self.param.dt;
    }

    fn setup(&mut self) {
        let cmaker = FccConfigMaker::new(0.5);
        cmaker.make_conf(&mut self.vars, &self.param);
        self.vars.set_initial_velocity(1.0);
        self.mesh.setup(&self.param, &self.vars);
    }

    pub fn run(&mut self) {
        self.setup();
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
