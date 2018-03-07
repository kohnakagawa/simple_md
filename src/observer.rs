use utils::f64vec3::F64vec3;

use parameter::Parameter;
use variables::Variables;
use meshlist::Pair;
use std::fs::File;
use std::io::BufWriter;

pub struct Observer {
    writer: BufWriter<File>,
}

impl Observer {
    pub fn new() -> Observer {
        Observer {
            writer: BufWriter::new(File::create("traject.xyz").unwrap()),
        }
    }

    pub fn traject(&mut self, vars: &Variables) {
        vars.export_xyz(&mut self.writer);
    }

    pub fn kinetic_energy(&mut self, vars: &Variables) -> f64 {
        vars.atoms.iter().fold(0.0, |k, atom|{
            k + atom.v.x * atom.v.x
              + atom.v.y * atom.v.y
              + atom.v.z * atom.v.z
        }) * 0.5 / (vars.number_of_atoms() as f64)
    }

    pub fn potential_energy(&mut self,
                            vars: &Variables,
                            pairs: &Vec<Pair>,
                            param: &Parameter) -> f64 {
        let cl2 = param.cl * param.cl;
        pairs.iter().fold(0.0, |u, pair| {
            let Pair(i, j) = *pair;
            let mut dr = F64vec3::new(
                vars.atoms[i].r.x - vars.atoms[j].r.x,
                vars.atoms[i].r.y - vars.atoms[j].r.y,
                vars.atoms[i].r.z - vars.atoms[j].r.z,
            );
            param.adjust_pbc(&mut dr);
            let r2 = dr.norm2();
            if r2 > cl2 {
                u
            } else {
                let r6 = r2 * r2 * r2;
                let r12 = r6 * r6;
                u + 4.0 * (1.0 / r12 - 1.0 / r6) + param.c0
            }
        }) / (vars.number_of_atoms() as f64)
    }
}
