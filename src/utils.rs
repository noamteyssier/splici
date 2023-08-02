use anyhow::{bail, Result};
use bedrs::{Coordinates, GenomicInterval};
use hashbrown::HashMap;
use noodles::core::{Position, Region};
use std::{hash::Hash, str::from_utf8};

/// Flips the K, V pairs in a HashMap
/// Used to reverse the name -> idx mappings to idx -> name
pub fn flip_map<K, V>(map: &HashMap<K, V>) -> HashMap<V, K>
where
    K: Eq + Hash + Clone,
    V: Eq + Hash + Clone,
{
    let mut flipped = HashMap::new();
    for (key, value) in map.iter() {
        flipped.insert(value.clone(), key.clone());
    }
    flipped
}

/// Converts a `GenomicInterval` to a `Region`
pub fn interval_to_region(
    giv: &GenomicInterval<usize>,
    genome_name_map: &HashMap<usize, Vec<u8>>,
) -> Result<Region> {
    let name = if let Some(name) = genome_name_map.get(&giv.chr()) {
        from_utf8(name)?.to_string()
    } else {
        bail!("Could not find genome name for id {}", giv.chr())
    };
    let start = Position::try_from(giv.start())?;
    let end = Position::try_from(giv.end())?;
    let region = Region::new(name, start..=end);
    Ok(region)
}
