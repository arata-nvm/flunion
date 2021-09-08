use std::collections::HashMap;

use crate::constants::*;
use crate::vecmath::vector2::Vector2;
use crate::Particle;

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct NeighborMapIndex(usize);
pub type NeighborMap = HashMap<NeighborMapIndex, Vec<usize>>;

pub fn new_neighbor_map(ps: &Vec<Particle>) -> NeighborMap {
    let mut map = NeighborMap::new();
    for (i, p) in ps.iter().enumerate() {
        insert_neighbor_map(&mut map, i, p.position);
    }
    map
}

fn insert_neighbor_map(map: &mut NeighborMap, i: usize, position: Vector2) {
    let index = calc_index(position);
    map.entry(index).or_insert_with(Vec::new).push(i);
}

pub fn neighbors(map: &NeighborMap, position: Vector2) -> Vec<usize> {
    let mut indices = Vec::new();
    let d = H / SPH_SIMSCALE;

    for x in -1..2 {
        for y in -1..2 {
            let v = Vector2::new(position[0] + x as f64 * d, position[1] + y as f64 * d);
            if (MIN[0] <= v[0] && v[0] <= MAX[0]) && (MIN[1] <= v[1] && v[1] <= MAX[1]) {
                let map_index = calc_index(v);
                if let Some(particles) = map.get(&map_index) {
                    indices.extend(particles);
                }
            }
        }
    }

    indices
}

const D: f64 = H / SPH_SIMSCALE;
const INV_D: f64 = 1.0 / D;
const MX: usize = ((MAX[0] - MIN[0]) * INV_D) as usize;

fn calc_index(position: Vector2) -> NeighborMapIndex {
    let x = ((position[0] - MIN[0]) * INV_D) as usize;
    let y = ((position[1] - MIN[1]) * INV_D) as usize;
    NeighborMapIndex(x + y * MX)
}
