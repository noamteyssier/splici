use std::str::Utf8Error;

use crate::{IdMap, NameMap};
use anyhow::Result;
use hashbrown::HashMap;

/// A struct to handle the translation of ids to names and vice versa
///
/// Used as an O(1) mapping between numeric ids and their corresponding
/// byte string names.
#[derive(Debug)]
pub struct Translater {
    /// A mapping from numeric ids to byte string names
    id_to_name: IdMap,

    /// A mapping from byte string names to numeric ids
    name_to_id: NameMap,
}
impl Translater {
    pub fn new() -> Self {
        Self {
            id_to_name: HashMap::new(),
            name_to_id: HashMap::new(),
        }
    }

    pub fn has_name(&self, name: &[u8]) -> bool {
        self.name_to_id.contains_key(name)
    }

    pub fn get_id_map(&self) -> &IdMap {
        &self.id_to_name
    }

    pub fn get_name_map(&self) -> &NameMap {
        &self.name_to_id
    }

    pub fn get_id_map_mut(&mut self) -> &mut IdMap {
        &mut self.id_to_name
    }

    pub fn get_name_map_mut(&mut self) -> &mut NameMap {
        &mut self.name_to_id
    }

    pub fn get_id(&self, name: &[u8]) -> Option<&usize> {
        self.name_to_id.get(name)
    }

    pub fn get_name(&self, id: usize) -> Option<&Vec<u8>> {
        self.id_to_name.get(&id)
    }

    pub fn get_name_str(&self, id: usize) -> Option<Result<&str, Utf8Error>> {
        self.id_to_name.get(&id).map(|v| std::str::from_utf8(v))
    }

    pub fn len(&self) -> usize {
        self.id_to_name.len()
    }

    pub fn insert(&mut self, name: Vec<u8>, id: usize) {
        self.id_to_name.insert(id, name.clone());
        self.name_to_id.insert(name, id);
    }
}
