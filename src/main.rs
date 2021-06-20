use flunion::vecmath::vector2::Vector2;
use flunion::{constants::*, step, Particle};

use tetra::{
    graphics::{
        self,
        mesh::{Mesh, ShapeStyle},
        Color, DrawParams,
    },
    math::Vec2,
    Context, ContextBuilder, State,
};

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
        step(&mut self.particles);

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
            let circle = Mesh::circle(ctx, ShapeStyle::Stroke(1.0), p, 20.0)?;
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
