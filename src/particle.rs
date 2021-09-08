use std::fs::File;
use std::io::Write;

use crate::vecmath::vector2;
use crate::vecmath::vector2::Vector2;
use crate::{constants::*, neighbors, new_neighbor_map, NeighborMap};

#[derive(Debug)]
pub struct Particle {
    pub position: Vector2,
    pub velocity: Vector2,
    pub force: Vector2,

    pub rho: f64,
    pub p: f64,
}

impl Particle {
    pub fn new(position: Vector2) -> Self {
        Self {
            position,
            velocity: Vector2::zero(),
            force: Vector2::zero(),

            rho: 0.0,
            p: 0.0,
        }
    }
}

pub fn step(ps: &mut Vec<Particle>) {
    let map = new_neighbor_map(ps);
    calculate_density_and_pressure(ps, &map);
    calculate_force(ps, &map);
    calculate_position(ps);
}

// 密度と圧力の計算
fn calculate_density_and_pressure(ps: &mut Vec<Particle>, map: &NeighborMap) {
    // 影響半径の距離
    let h2 = H * H;

    for i in 0..ps.len() {
        let mut sum = 0.0;
        let p1 = &ps[i];
        let neighbors = neighbors(map, p1.position);
        for j in neighbors {
            if i == j {
                continue;
            }

            let p2 = &ps[j];

            // 2粒子間の距離
            let dr = p1.position.sub(p2.position).mul(SPH_SIMSCALE);
            let r2 = dr.len_sq();

            // 粒子が影響半径内にある
            if h2 > r2 {
                let c = h2 - r2;
                sum += c * c * c;
            }
        }

        let mut p = &mut ps[i];
        p.rho = sum * SPH_PMASS * *POLY6KERN; // 密度
        p.p = (p.rho - SPH_RESTDENSITY) * SPH_INTSTIFF; // 圧力
        p.rho = 1.0 / p.rho;
    }
}

// 力の計算
fn calculate_force(ps: &mut Vec<Particle>, map: &NeighborMap) {
    for i in 0..ps.len() {
        let mut force = Vector2::zero();
        let p1 = &ps[i];
        let neighbors = neighbors(map, p1.position);
        for j in neighbors {
            if i == j {
                continue;
            }

            let p2 = &ps[j];

            // 2粒子間の距離
            let dr = p1.position.sub(p2.position).mul(SPH_SIMSCALE);
            let r = dr.len();

            // 粒子が影響半径内にある
            if H > r {
                let c = H - r;
                let pterm = -0.5 * c * *SPIKYKERN * (p1.p + p2.p) / r; // 圧力項
                let vterm = *LAPKERN * SPH_VISC; // 粘性項
                let mut fcurr = dr.mul(pterm).add(p2.velocity.sub(p1.velocity).mul(vterm));
                fcurr.set_mul(c * p1.rho * p2.rho);
                force.set_add(fcurr);
            }
        }

        ps[i].force = force;
    }
}

// 位置更新
// 境界条件を適用
fn calculate_position(ps: &mut Vec<Particle>) {
    let g = Vector2::new(0.0, -9.8);

    for p in ps {
        let mut accel = p.force.mul(SPH_PMASS);

        let speed = accel.len_sq();
        if speed > SPH_LIMIT * SPH_LIMIT {
            accel.set_mul(SPH_LIMIT / speed.sqrt());
        }

        // --- X ---

        let diff = 2.0 * SPH_RADIUS - (p.position[0] - MIN[0]) * SPH_SIMSCALE;
        if diff > SPH_EPSILON {
            let norm = Vector2::new(1.0, 0.0);
            let adj = SPH_EXTSTIFF * diff - SPH_EXTDAMP * vector2::dot(norm, p.velocity);
            accel.set_add(norm.mul(adj));
        }

        let diff = 2.0 * SPH_RADIUS - (MAX[0] - p.position[0]) * SPH_SIMSCALE;
        if diff > SPH_EPSILON {
            let norm = Vector2::new(-1.0, 0.0);
            let adj = SPH_EXTSTIFF * diff - SPH_EXTDAMP * vector2::dot(norm, p.velocity);
            accel.set_add(norm.mul(adj));
        }

        // --- X ---

        // --- Y ---

        let diff = 2.0 * SPH_RADIUS - (p.position[1] - MIN[1]) * SPH_SIMSCALE;
        if diff > SPH_EPSILON {
            let norm = Vector2::new(0.0, 1.0);
            let adj = SPH_EXTSTIFF * diff - SPH_EXTDAMP * vector2::dot(norm, p.velocity);
            accel.set_add(norm.mul(adj));
        }

        let diff = 2.0 * SPH_RADIUS - (MAX[1] - p.position[1]) * SPH_SIMSCALE;
        if diff > SPH_EPSILON {
            let norm = Vector2::new(0.0, -1.0);
            let adj = SPH_EXTSTIFF * diff - SPH_EXTDAMP * vector2::dot(norm, p.velocity);
            accel.set_add(norm.mul(adj));
        }

        // --- Y ---

        accel.set_add(g);
        p.velocity.set_add(accel.mul(DT));
        p.position.set_add(p.velocity.mul(DT).div(SPH_SIMSCALE));
    }
}

pub fn output_particles(i: usize, ps: &Vec<Particle>) {
    let file_name = format!("result{:03}.pov", i);
    println!("processing {} ...", file_name);

    let mut s = String::from(
        "#include \"colors.inc\"\ncamera { location <10, 30, -40.0> look_at <10, 10, 0.0> }\nlight_source { <0, 30, -30> color White }\n",
    );

    for p in ps {
        s.push_str("sphere {\n");
        s.push_str(&format!(
            "  <{}, {}, 0>, 0.5\n",
            p.position[0], p.position[1]
        ));
        s.push_str("  texture { pigment { color Gray30 } }\n}\n");
    }

    let mut f = File::create(file_name).unwrap();
    write!(f, "{}", s).unwrap();
}
