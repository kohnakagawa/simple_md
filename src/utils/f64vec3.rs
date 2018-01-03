#[derive(Debug, Clone, Copy)]
pub struct F64vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl F64vec3 {
    pub fn new(x: f64, y: f64, z: f64) -> F64vec3 {
        F64vec3 {x: x, y: y, z: z}
    }

    pub fn zero() -> F64vec3 {
        F64vec3 {x: 0.0, y: 0.0, z: 0.0}
    }

    pub fn norm2(&self) -> f64 {
        self.x*self.x + self.y*self.y + self.z*self.z
    }
}
