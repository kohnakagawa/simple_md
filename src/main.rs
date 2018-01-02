// A simple O(N^2) molecular dynamics simulation code for my study.
// step1
extern crate rand;

use std::f64;
use std::io::{BufWriter, Write};
use rand::distributions::{IndependentSample, Range};

#[derive(Debug, Clone, Copy)]
struct F64vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl F64vec3 {
    fn new(x: f64, y: f64, z: f64) -> F64vec3 {
        F64vec3 {x: x, y: y, z: z}
    }

    fn zero() -> F64vec3 {
        F64vec3 {x: 0.0, y: 0.0, z: 0.0}
    }

    fn norm2(&self) -> f64 {
        self.x*self.x + self.y*self.y + self.z*self.z
    }
}

#[derive(Debug, Clone, Copy)]
struct Atom {
    r: F64vec3,
    v: F64vec3,
}

impl Atom {
    fn new (r: F64vec3, v: F64vec3) -> Atom {
        Atom {r: r, v: v}
    }
}

struct Variables {
    atoms: Vec<Atom>,
    time: f64,
}

impl Variables {
    fn new() -> Variables {
        Variables {atoms: Vec::new(),
                   time: 0.0}
    }

    fn add_atoms(&mut self, r: F64vec3) {
        self.atoms.push(Atom::new(r, F64vec3::zero()));
    }

    fn export_xyz<W: Write>(&self, f: &mut BufWriter<W>) {
        f.write_fmt(format_args!("{}\n", self.atoms.len()))
            .unwrap();
        f.write(b"xyz file format\n").unwrap();
        for atom in &self.atoms {
            f.write_fmt(format_args!("O {:.10} {:.10} {:.10}\n",
                                     atom.r.x, atom.r.y, atom.r.z))
                .unwrap();
        }
    }

    fn number_of_atoms(&self) -> usize {
        self.atoms.len()
    }

    fn set_initial_velocity(&mut self, v0: f64) {
        let mut rng = rand::thread_rng();
        let between = Range::new(0.0, 1.0);
        let mut av  = F64vec3::new(0.0, 0.0, 0.0);
        for atom in &mut self.atoms {
            let z = between.ind_sample(&mut rng);
            let phi = 2.0 * between.ind_sample(&mut rng) * f64::consts::PI;
            let v_ini =
                F64vec3::new(v0 * (1.0 - z * z).sqrt() * phi.cos(),
                             v0 * (1.0 - z * z).sqrt() * phi.sin(),
                             v0 * z);
            atom.v = v_ini;
            av.x += v_ini.x;
            av.y += v_ini.y;
            av.z += v_ini.z;
        }
        let inv_pn = 1.0 / (self.atoms.len() as f64);
        av.x *= inv_pn;
        av.y *= inv_pn;
        av.z *= inv_pn;
        for atom in &mut self.atoms {
            atom.v.x -= av.x;
            atom.v.y -= av.y;
            atom.v.z -= av.z;
        }
    }
}

#[derive(Debug)]
struct Parameter {
    box_len: F64vec3,
    dt: f64,
    cl: f64,
    c0: f64,
}

impl Parameter {
    fn new() -> Parameter {
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
        }
    }

    fn adjust_pbc(&self, dr: &mut F64vec3) {
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
        if dr.z <  bh.z {dr.z -= bl.z;}
    }
}

struct Observer {}

impl Observer {
    fn new() -> Observer {
        Observer {}
    }

    fn kinetic_energy(&mut self, vars: &Variables) -> f64 {
        let mut k = 0.0;
        for atom in &vars.atoms {
            k += atom.v.x * atom.v.x;
            k += atom.v.y * atom.v.y;
            k += atom.v.z * atom.v.z;
        }
        k /= vars.number_of_atoms() as f64;
        k * 0.5
    }

    fn potential_energy(&mut self, vars: &Variables, param: &Parameter) -> f64 {
        let mut u = 0.0;
        let n_atoms = vars.number_of_atoms();
        let cl2 = param.cl * param.cl;
        for i in 0..n_atoms-1 {
            for j in i+1..n_atoms {
                let mut dr = F64vec3::new(
                    vars.atoms[i].r.x - vars.atoms[j].r.x,
                    vars.atoms[i].r.y - vars.atoms[j].r.y,
                    vars.atoms[i].r.z - vars.atoms[j].r.z,
                );
                param.adjust_pbc(&mut dr);
                let r2 = dr.norm2();
                if r2 > cl2 {continue;}
                let r6 = r2 * r2 * r2;
                let r12 = r6 * r6;
                u += 4.0 * (1.0 / r12 - 1.0 / r6) + param.c0;
            }
        }
        u /= n_atoms as f64;
        u
    }
}

struct ConfigMaker {}

impl ConfigMaker {
    fn make_fcc(vars: &mut Variables, param: &Parameter) {
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

struct MDSystem {
    vars: Variables,
    param: Parameter,
    obs: Observer,
}

impl MDSystem {
    fn new() -> MDSystem {
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
        for i in 0..n_atoms-2 {
            for j in i+1..n_atoms-1 {
                let mut dr = F64vec3::new(
                    vars.atoms[i].r.x - vars.atoms[j].r.x,
                    vars.atoms[i].r.y - vars.atoms[j].r.y,
                    vars.atoms[i].r.z - vars.atoms[j].r.z,
                );
                param.adjust_pbc(&mut dr);
                let r2 = dr.norm2();
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

    fn run(&mut self) {
        ConfigMaker::make_fcc(&mut self.vars,
                              &self.param);
        self.vars.set_initial_velocity(1.0);
        let steps = 10000;
        let observe = 100;
        for i in 0..steps {
            if (i % observe) == 0 {
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

fn main() {
    let mut mdsystem = MDSystem::new();
    mdsystem.run();
}

