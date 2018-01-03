extern crate rand;

use utils::f64vec3::F64vec3;
use self::rand::distributions::{IndependentSample, Range};
use std::io::{BufWriter, Write};
use std::f64;

#[derive(Debug, Clone, Copy)]
pub struct Atom {
    pub r: F64vec3,
    pub v: F64vec3,
}

impl Atom {
    pub fn new (r: F64vec3, v: F64vec3) -> Atom {
        Atom {r: r, v: v}
    }
}

pub struct Variables {
    pub atoms: Vec<Atom>,
    pub time: f64,
}

impl Variables {
    pub fn new() -> Variables {
        Variables {atoms: Vec::new(),
                   time: 0.0}
    }

    pub fn add_atoms(&mut self, r: F64vec3) {
        self.atoms.push(Atom::new(r, F64vec3::zero()));
    }

    pub fn export_xyz<W: Write>(&self, f: &mut BufWriter<W>) {
        f.write_fmt(format_args!("{}\n", self.atoms.len()))
            .unwrap();
        f.write(b"xyz file format\n").unwrap();
        for atom in &self.atoms {
            f.write_fmt(format_args!("O {:.10} {:.10} {:.10} {:.10} {:.10} {:.10}\n",
                                     atom.r.x, atom.r.y, atom.r.z,
                                     atom.v.x, atom.v.y, atom.v.z))
                .unwrap();
        }
    }

    pub fn number_of_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn set_initial_velocity(&mut self, v0: f64) {
        let mut rng = rand::thread_rng();
        let between = Range::new(0.0, 1.0);
        let mut av  = F64vec3::new(0.0, 0.0, 0.0);
        for atom in &mut self.atoms {
            let z = between.ind_sample(&mut rng) * 2.0 - 1.0;
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
