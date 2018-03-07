use utils::f64vec3::F64vec3;
use std::path::Path;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug)]
pub struct Parameter {
    pub box_len: F64vec3,
    pub dt: f64,
    pub cl: f64,
    pub c0: f64,
    pub margin: f64,
    pub steps: u32,
    pub observe: u32,
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
            steps: 100000,
            observe: 1000,
        }
    }

    pub fn read_from_file(&mut self, fname: &Path) {
        let reader = BufReader::new(File::open(fname).unwrap());
        let lines = reader.lines();
        for line in lines {
            let line = line.unwrap();
            let mut words: Vec<&str> = line.split_whitespace().collect();
            if words[0] == "cl" {self.cl = words[1].parse().unwrap();}
            if words[0] == "dt" {self.dt = words[1].parse().unwrap();}
            if words[0] == "margin" {self.margin = words[1].parse().unwrap();}
            if words[0] == "steps" {self.steps  = words[1].parse().unwrap();}
            if words[0] == "observe" {self.observe = words[1].parse().unwrap();}
            if words[0] == "box_len.x" {self.box_len.x = words[1].parse().unwrap();}
            if words[0] == "box_len.y" {self.box_len.y = words[1].parse().unwrap();}
            if words[0] == "box_len.z" {self.box_len.z = words[1].parse().unwrap();}
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
