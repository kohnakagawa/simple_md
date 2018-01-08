use variables::Variables;
use parameter::Parameter;
use meshlist::MeshList;
use meshlist::Pair;
use utils::f64vec3::F64vec3;

pub struct ForceCalculator {}

impl ForceCalculator {
    pub fn calculate_force(vars: &mut Variables,
                           meshlist: &MeshList,
                           param: &Parameter) {
        // ForceCalculator::
        // calculate_brute_force(vars, param);

        // ForceCalculator::
        // calculate_pair(vars, &meshlist.pairs, param);

        ForceCalculator::
        calculate_sorted_list(vars, &meshlist.pointer,
                              &meshlist.sorted_list,
                              param);

    }

    #[allow(dead_code)]
    fn calculate_brute_force(vars: &mut Variables,
                             param: & Parameter) {
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

    #[allow(dead_code)]
    fn calculate_pair(vars: &mut Variables,
                      pairs: &Vec<Pair>,
                      param: &Parameter) {
        let cl2 = param.cl * param.cl;
        let dt  = param.dt;
        for p in pairs {
            let Pair(i, j) = *p;
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

    #[allow(dead_code)]
    fn calculate_sorted_list(vars: &mut Variables,
                             pointer: &Vec<usize>,
                             sorted_list: &Vec<usize>,
                             param: &Parameter) {
        let n_atoms = vars.number_of_atoms();
        let cl2 = param.cl * param.cl;
        let dt  = param.dt;
        for i in 0..n_atoms {
            let ri = vars.atoms[i].r;
            let mut pi = F64vec3::zero();
            for k in pointer[i]..pointer[i+1] {
                let j = sorted_list[k];
                let mut dr = F64vec3::new(
                    vars.atoms[j].r.x - ri.x,
                    vars.atoms[j].r.y - ri.y,
                    vars.atoms[j].r.z - ri.z,
                );
                param.adjust_pbc(&mut dr);
                let r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
                if r2 > cl2 {continue;}
                let r6 = r2 * r2 * r2;
                let df = (24.0 * r6 - 48.0) / (r6 * r6 * r2) * dt;
                pi.x += df * dr.x;
                pi.y += df * dr.y;
                pi.z += df * dr.z;
                vars.atoms[j].v.x -= df * dr.x;
                vars.atoms[j].v.y -= df * dr.y;
                vars.atoms[j].v.z -= df * dr.z;
            }
            vars.atoms[i].v.x += pi.x;
            vars.atoms[i].v.y += pi.y;
            vars.atoms[i].v.z += pi.z;
        }

    }
}
