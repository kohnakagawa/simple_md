// A simple O(N^2) molecular dynamics simulation code for my study.
// step1
extern crate rand;

use std::f64;
use std::fs::File;
use std::io::{BufWriter, Write};
use rand::distributions::{IndependentSample, Range};

#[derive(Debug, Clone, Copy)]
struct F64vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl F64vec3 {
    fn new(x_: f64, y_: f64, z_: f64) -> F64vec3 {
        F64vec3 {x: x_, y: y_, z: z_}
    }

    fn zeros() -> F64vec3 {
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
    fn new (r_: F64vec3, v_: F64vec3) -> Atom {
        Atom {r: r_, v: v_}
    }
}

struct Variables {
    atoms: Vec<Atom>,
    time: f64,
}

impl Variables {
    fn add_atoms(&mut self, r: F64vec3) {
        self.atoms.push(Atom::new(r, F64vec3::zeros()));
    }

    fn export_xyz<W: Write>(&self, f: &mut BufWriter<W>) {
        f.write_fmt(format_args!("{}\n", self.atoms.len()))
            .unwrap();
        f.write("xyz file format\n").unwrap();
        for atom in self.atoms {
            f.write_fmt(format_args!("O {:.10} {:.10} {:.10}\n",
                                     atom.r.x, atom.r.y, atom.r.z))
                .unwrap();
        }
    }

    fn number_of_atoms(&self) -> usize {
        self.atoms.len()
    }

    fn set_initial_velocity(&self, v0: f64) {
        let mut rng = rand::thread_rng();
        let between = Range::new(0.0, 1.0);
        let av      = F64vec3::new(0.0, 0.0, 0.0);
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
    fn kinetic_energy(vars: &Variables) -> f64 {
        let mut k = 0.0;
        for atom in &vars.atoms {
            k += atom.v.x * atom.v.x;
            k += atom.v.y * atom.v.y;
            k += atom.v.z * atom.v.z;
        }
        k /= vars.number_of_atoms() as f64;
        k * 0.5
    }

    fn potential_energy(vars: &Variables, param: &Parameter) -> f64 {
        let mut u = 0.0;
        let n_atoms = vars.number_of_atoms();
        let cl2 = param.cl * param.cl;
        for i in 0..n_atoms-2 {
            for j in i+1..n_atoms-1 {
                let dr = F64vec3::new(
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
    fn make_fcc() {
        let density = 0.5;

    }
}

fn main() {

}
