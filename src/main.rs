use flunion::vecmath::{
    self,
    vector2::{self, Vector2},
};
use once_cell::unsync::Lazy;
use std::{
    f64::consts::PI,
    fs::File,
    io::Write,
    ops::{Mul, Sub},
};
use tetra::{
    graphics::{
        self,
        mesh::{Mesh, ShapeStyle},
        Color, DrawParams,
    },
    math::Vec2,
    Context, ContextBuilder, State,
};

const MAX_LOOP: usize = 1000;
const H: f64 = 0.01;
const SPH_PMASS: f64 = 0.00020543;
const SPH_INTSTIFF: f64 = 1.00;
const SPH_EXTSTIFF: f64 = 10000.0;
const SPH_EXTDAMP: f64 = 256.0;
const POLY6KERN: Lazy<f64> = Lazy::new(|| 315.0 / (64.0 * PI * H.powf(9.0)));
const SPIKYKERN: Lazy<f64> = Lazy::new(|| -45.0 / (PI * H.powf(6.0)));
const LAPKERN: Lazy<f64> = Lazy::new(|| 45.0 / (PI * H.powf(6.0)));
const DT: f64 = 0.004;
const RADIUS: f64 = 0.004;
const EPSILON: f64 = 0.00001;
const INITMIN: [f64; 3] = [0.0, 0.0, 0.0];
const INITMAX: [f64; 3] = [10.0, 20.0, 0.0];
const MIN: [f64; 3] = [0.0, 0.0, -10.0];
const MAX: [f64; 3] = [20.0, 20.0, 10.0];
const SPH_RESTDENSITY: f64 = 600.0;
const SPH_PDIST: Lazy<f64> = Lazy::new(|| (SPH_PMASS / SPH_RESTDENSITY).powf(1.0 / 3.0));
const SPH_SIMSCALE: f64 = 0.004;
const SPH_VISC: f64 = 0.2;
const SPH_LIMIT: f64 = 200.0;
const D: Lazy<f64> = Lazy::new(|| *SPH_PDIST * 0.87 / SPH_SIMSCALE);

#[derive(Debug)]
struct Particle {
    position: Vector2,
    velocity: Vector2,
    force: Vector2,

    rho: f64,
    p: f64,
}

impl Particle {
    fn new(position: Vector2) -> Self {
        Self {
            position,
            velocity: Vector2::zero(),
            force: Vector2::zero(),

            rho: 0.0,
            p: 0.0,
        }
    }
}

// 密度と圧力の計算
fn calculate_density_and_pressure(ps: &mut Vec<Particle>) {
    // 影響半径の距離
    let h2 = H * H;

    for i in 0..ps.len() {
        let mut sum = 0.0;
        for j in 0..ps.len() {
            if i == j {
                continue;
            }

            let p1 = &ps[i];
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
fn calculate_force(ps: &mut Vec<Particle>) {
    for i in 0..ps.len() {
        let mut force = Vector2::zero();
        for j in 0..ps.len() {
            if i == j {
                continue;
            }

            let p1 = &ps[i];
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

        let diff = 2.0 * RADIUS - (p.position[0] - MIN[0]) * SPH_SIMSCALE;
        if diff > EPSILON {
            let norm = Vector2::new(1.0, 0.0);
            let adj = SPH_EXTSTIFF * diff - SPH_EXTDAMP * vector2::dot(norm, p.velocity);
            accel.set_add(norm.mul(adj));
        }

        let diff = 2.0 * RADIUS - (MAX[0] - p.position[0]) * SPH_SIMSCALE;
        if diff > EPSILON {
            let norm = Vector2::new(-1.0, 0.0);
            let adj = SPH_EXTSTIFF * diff - SPH_EXTDAMP * vector2::dot(norm, p.velocity);
            accel.set_add(norm.mul(adj));
        }

        // --- X ---

        // --- Y ---

        let diff = 2.0 * RADIUS - (p.position[1] - MIN[1]) * SPH_SIMSCALE;
        if diff > EPSILON {
            let norm = Vector2::new(0.0, 1.0);
            let adj = SPH_EXTSTIFF * diff - SPH_EXTDAMP * vector2::dot(norm, p.velocity);
            accel.set_add(norm.mul(adj));
        }

        let diff = 2.0 * RADIUS - (MAX[1] - p.position[1]) * SPH_SIMSCALE;
        if diff > EPSILON {
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

fn output_particles(i: usize, ps: &Vec<Particle>) {
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

fn init() -> Vec<Particle> {
    let mut v = Vec::new();

    let mut y = INITMIN[1];
    while y <= INITMAX[1] {
        let mut x = INITMIN[0];
        while x <= INITMAX[0] {
            println!("{}, {}", x, y);
            v.push(Particle::new(Vector2::new(x, y)));
            x += *D;
        }
        y += *D;
    }

    v
}

struct GameState {
    particles: Vec<Particle>,
}

impl GameState {
    fn new(_ctx: &mut Context) -> tetra::Result<Self> {
        Ok(Self { particles: init() })
    }
}

impl State for GameState {
    fn update(&mut self, ctx: &mut Context) -> tetra::Result {
        calculate_density_and_pressure(&mut self.particles);
        calculate_force(&mut self.particles);
        calculate_position(&mut self.particles);

        tetra::window::set_title(ctx, format!("fluid fps:{:.02}", tetra::time::get_fps(ctx)));
        Ok(())
    }

    fn draw(&mut self, ctx: &mut Context) -> tetra::Result {
        graphics::clear(ctx, Color::BLACK);

        for p in &self.particles {
            let scale = 800.0 / 20.0;
            let p = Vec2::new(
                p.position[0] as f32 * scale,
                800.0 - p.position[1] as f32 * scale,
            );
            let circle = Mesh::circle(ctx, ShapeStyle::Stroke(1.0), p, 10.0)?;
            circle.draw(ctx, DrawParams::new());
        }

        Ok(())
    }
}

fn main() -> tetra::Result {
    ContextBuilder::new("fluid", 800, 800)
        .build()?
        .run(GameState::new)
}
