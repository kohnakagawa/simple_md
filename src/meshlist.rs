use variables::Variables;
use parameter::Parameter;
use utils::f64vec3::F64vec3;

const NUM_NEIGHBOR_MESH: usize = 13;

pub struct Pair(pub usize, pub usize);

pub struct MeshList {
    mesh_size: F64vec3,
    imesh_size: F64vec3,
    sl: f64,
    mx: i32,
    my: i32,
    mz: i32,
    number_of_mesh: usize,
    neighbor_mesh_id: Vec<[usize; NUM_NEIGHBOR_MESH]>,
    count: Vec<usize>,
    indexes: Vec<usize>,
    sorted_buffer: Vec<usize>,
    particle_position: Vec<usize>,
    ptcl_id_of_nmesh: Vec<Vec<usize>>,
    num_neighbor: Vec<usize>,
    pub pairs: Vec<Pair>,
    pub pointer: Vec<usize>,
    pointer2: Vec<usize>,
    pub sorted_list: Vec<usize>,
}

impl MeshList {
    pub fn new() -> MeshList {
        MeshList {
            mesh_size: F64vec3::zero(),
            imesh_size: F64vec3::zero(),
            sl: 0.0,
            mx: 0,
            my: 0,
            mz: 0,
            number_of_mesh: 0,
            neighbor_mesh_id: Vec::new(),
            count: Vec::new(),
            indexes: Vec::new(),
            sorted_buffer: Vec::new(),
            particle_position: Vec::new(),
            ptcl_id_of_nmesh: Vec::new(),
            num_neighbor: Vec::new(),
            pairs: Vec::new(),
            pointer: Vec::new(),
            pointer2: Vec::new(),
            sorted_list: Vec::new(),
        }
    }

    pub fn setup(&mut self, param: &Parameter, vars: &Variables) {
        let sl = param.cl + param.margin;
        let (mx, my, mz) = (
            (param.box_len.x / sl) as i32,
            (param.box_len.y / sl) as i32,
            (param.box_len.z / sl) as i32,
        );
        let (mlx, mly, mlz) = (
            param.box_len.x / mx as f64,
            param.box_len.y / my as f64,
            param.box_len.z / mz as f64,
        );

        assert!(mx > 1, "mx should be > 1");
        assert!(my > 1, "my should be > 1");
        assert!(mz > 1, "mz should be > 1");

        self.mesh_size = F64vec3::new(mlx, mly, mlz);
        self.imesh_size = F64vec3::new(1.0 / mlx, 1.0 / mly, 1.0 / mlz);
        self.sl = sl;
        self.mx = mx; self.my = my; self.mz = mz;
        self.number_of_mesh    = (mx*my*mz) as usize;
        self.neighbor_mesh_id  = vec![[0;NUM_NEIGHBOR_MESH];self.number_of_mesh];
        self.count             = vec![0; self.number_of_mesh];
        self.indexes           = vec![0; self.number_of_mesh + 1];
        self.sorted_buffer     = vec![0; vars.number_of_atoms()];
        self.particle_position = vec![0; vars.number_of_atoms()];
        self.ptcl_id_of_nmesh  = vec![vec![]; self.number_of_mesh];
        self.num_neighbor      = vec![0; vars.number_of_atoms()];
        self.pointer           = vec![0; vars.number_of_atoms() + 1];
        self.pointer2          = vec![0; vars.number_of_atoms() + 1];

        self.set_neighbor_mesh_id();
    }

    fn apply_pbc(&self, idx: &mut (i32, i32, i32)) {
        if idx.0 <  0       {idx.0 += self.mx}
        if idx.0 >= self.mx {idx.0 -= self.mx}
        if idx.1 <  0       {idx.1 += self.my}
        if idx.1 >= self.my {idx.1 -= self.my}
        if idx.2 <  0       {idx.2 += self.mz}
        if idx.2 >= self.mz {idx.2 -= self.mz}
    }

    fn set_neighbor_mesh_id_sub(&mut self, ix: i32, iy: i32, iz: i32, icnt: usize) {
        let mut jcnt = 0;
        for jz in -1..2 { for jy in -1..2 { for jx in -1..2 {
            let mut idx = (
                ix + jx,
                iy + jy,
                iz + jz,
            );
            self.apply_pbc(&mut idx);
            self.neighbor_mesh_id[icnt][jcnt]
                = self.mesh_id(idx.0, idx.1, idx.2);
            jcnt += 1;
            if jcnt == NUM_NEIGHBOR_MESH {return;}
        }}}
    }

    fn set_neighbor_mesh_id(&mut self) {
        let mut icnt = 0;
        for iz in 0..self.mz { for iy in 0..self.my { for ix in 0..self.mx {
            self.set_neighbor_mesh_id_sub(ix, iy, iz, icnt);
            icnt += 1;
        }}}
    }

    fn make_ptcl_id_of_nmesh(&mut self) {
        for pid in &mut self.ptcl_id_of_nmesh { pid.clear(); }
        for imesh in 0..self.number_of_mesh {
            let beg = self.indexes[imesh    ];
            let end = self.indexes[imesh + 1];
            self.ptcl_id_of_nmesh[imesh]
                .extend_from_slice(&self.sorted_buffer[beg..end]);
            for jmesh in 0..NUM_NEIGHBOR_MESH {
                let nmesh_id = self.neighbor_mesh_id[imesh][jmesh];
                let beg = self.indexes[nmesh_id    ];
                let end = self.indexes[nmesh_id + 1];
                self.ptcl_id_of_nmesh[imesh]
                    .extend_from_slice(&self.sorted_buffer[beg..end]);
            }
        }
    }

    fn register_pair(&mut self, i: usize, j: usize) {
        let key: usize;
        let partner: usize;
        if i < j {
            key = i;
            partner = j;
        } else {
            key = j;
            partner = i;
        }
        self.pairs.push(Pair(key, partner));
    }

    fn search_mesh(&mut self, vars: &Variables, param: &Parameter) {
        self.pairs.clear();
        let sl2 = self.sl * self.sl;
        self.make_ptcl_id_of_nmesh();
        for imesh in 0..self.number_of_mesh {
            let ibeg = self.indexes[imesh];
            let icnt = self.count[imesh];
            for i in 0..icnt {
                let i_id = self.sorted_buffer[i + ibeg];
                let j_num = self.ptcl_id_of_nmesh[imesh].len();
                for j in i+1..j_num {
                    let j_id = self.ptcl_id_of_nmesh[imesh][j];
                    let ri = vars.atoms[i_id].r;
                    let rj = vars.atoms[j_id].r;
                    let mut dr = F64vec3::new(
                        ri.x - rj.x,
                        ri.y - rj.y,
                        ri.z - rj.z,
                    );
                    param.adjust_pbc(&mut dr);
                    let r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
                    if r2 < sl2 {self.register_pair(i_id, j_id)}
                }
            }
        }
    }

    // for debug
    #[allow(dead_code)]
    fn search_brute_force(&mut self, vars: &Variables, param: &Parameter) {
        self.pairs.clear();
        let sl2 = self.sl * self.sl;
        let num_atoms = vars.number_of_atoms();
        for i_id in 0..num_atoms-1 {
            for j_id in i_id+1..num_atoms {
                let mut dr = F64vec3::new(
                    vars.atoms[j_id].r.x - vars.atoms[i_id].r.x,
                    vars.atoms[j_id].r.y - vars.atoms[i_id].r.y,
                    vars.atoms[j_id].r.z - vars.atoms[i_id].r.z,
                );
                param.adjust_pbc(&mut dr);
                let r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
                if r2 < sl2 {self.register_pair(i_id, j_id)}
            }
        }

        println!("{} pairs are registered", self.pairs.len());
    }

    fn mesh_id(&self, ix: i32, iy: i32, iz: i32) -> usize {
        (ix + self.mx * (iy + iz * self.my)) as usize
    }

    fn get_mesh_pos(&self, r: &F64vec3) -> usize {
        let mut idx = (
            (r.x * self.imesh_size.x) as i32,
            (r.y * self.imesh_size.y) as i32,
            (r.z * self.imesh_size.z) as i32,
        );
        self.apply_pbc(&mut idx);
        self.mesh_id(idx.0, idx.1, idx.2)
    }

    fn register_mesh_pos(&mut self, vars: &Variables) {
        self.count = self.count.iter().map(|_| 0).collect();
        for (i, atom) in vars.atoms.iter().enumerate() {
            let mesh_pos = self.get_mesh_pos(&atom.r);
            self.particle_position[i] = mesh_pos;
            self.count[mesh_pos] += 1;
        }
        self.indexes[0] = 0;
        for i in 0..self.number_of_mesh {
            self.indexes[i+1] = self.indexes[i] + self.count[i];
        }
        assert!(self.indexes[self.number_of_mesh] == vars.number_of_atoms(),
                "prefix_sum is not correctly performed.");

        self.count = self.count.iter().map(|_| 0).collect();
        for i in 0..vars.number_of_atoms() {
            let mesh_pos                  = self.particle_position[i];
            let cnt                       = self.count[mesh_pos];
            let beg                       = self.indexes[mesh_pos];
            self.sorted_buffer[cnt + beg] = i;
            self.count[mesh_pos]          = cnt + 1;
        }

        self.check_sorted(&vars);
    }

    fn check_sorted(&self, vars: &Variables) {
        for i in 0..self.number_of_mesh {
            let beg = self.indexes[i];
            let end = beg + self.count[i];
            for k in beg..end {
                let j = self.sorted_buffer[k];
                assert!(self.get_mesh_pos(&vars.atoms[j].r) == i,
                        "mesh position is correctly registered.");
            }
        }
    }

    fn make_sorted_list_from_pairs(&mut self, vars: &Variables) {
        if self.sorted_list.len() < self.pairs.len() {
            self.sorted_list.resize(self.pairs.len() * 2, 0);
        }

        self.num_neighbor = self.num_neighbor.iter().map(|_| 0).collect();
        for p in &self.pairs {
            self.num_neighbor[p.0] += 1;
        }
        self.pointer[0] = 0;
        self.pointer2[0] = 0;
        for i in 0..vars.number_of_atoms() {
            self.pointer[i+1]  = self.pointer[i] + self.num_neighbor[i];
            self.pointer2[i+1] = self.pointer[i] + self.num_neighbor[i];
        }

        debug_assert!(self.pointer[vars.number_of_atoms()] == self.pairs.len(),
                      "prefix sum is not correctly performed.");

        for p in &self.pairs {
            let Pair(key, partner) = *p;
            let dst = self.pointer2[key];
            self.sorted_list[dst] = partner;
            self.pointer2[key] = dst + 1;
        }
    }

    pub fn make_sorted_list(&mut self, vars: &Variables, param: &Parameter) {
        self.register_mesh_pos(vars);
        self.make_ptcl_id_of_nmesh();
        self.search_mesh(vars, param);
        self.make_sorted_list_from_pairs(vars);
    }
}
