use std::fs::File;
use std::io::{self, Read, Write};
use nalgebra as na;
use fnv::FnvBuildHasher;
use by_addr::*;

fn split2_inclusive(line: &str, pat: &str) -> Vec<String> {
    match line.find(pat) {
        Some(entry) => {
            let first = line[..=entry].trim().to_owned();
            let last_str = line[entry + 1 ..].trim();
            if last_str.is_empty() {
                vec![first]
            } else {
                vec![first, last_str.to_owned()]
            }
        },
        None => vec![line.to_owned()],
    }    
}

fn str_between<'a>(line: &'a str, start_pat: &str, end_pat: &str) -> Option<&'a str> {
    let start = line.find(start_pat)? + start_pat.len();
    let end = line.find(end_pat)?;
    Some(&line[start..end])
}

#[derive(Debug)]
pub struct GeoFile(File);

impl GeoFile {
    pub fn open(path: &str) -> io::Result<Self> {
        File::open(path).map(GeoFile)
    }

    pub fn read_exprs(&mut self) -> Vec<String> {
        let mut buf = String::new();
        self.0.read_to_string(&mut buf).unwrap();
        buf.lines().flat_map(|x| split2_inclusive(x, ";")).map(|x| x.trim().to_owned()).collect()
    }
}

#[derive(Debug)]
pub struct OccFile(File);

impl OccFile {
    pub fn create(path: &str) -> io::Result<Self> {
        let mut file = File::create(path)?;
        writeln!(&mut file, "SetFactory(\"OpenCASCADE\");\nGeometry.OCCAutoFix = 0;").unwrap();
        Ok(OccFile(file))
    }

    fn write_point(&mut self, tag: u64, pos: &na::Point3<f64>) -> io::Result<()> {
        writeln!(self.0, "Point({})={{{},{},{}}};", tag, pos[0], pos[1], pos[2])?;
        Ok(())
    }

    fn write_line(&mut self, tag: u64, points: &[u64; 2]) -> io::Result<()> {
        writeln!(self.0, "Line({})={{{},{}}};", tag, points[0], points[1])?;
        Ok(())
    }

    fn write_line_loop(&mut self, tag: u64, lines: &[i64]) -> io::Result<()> {
        writeln!(self.0, "Line Loop({})={{{}}};", 
            tag, 
            lines.iter().map(|&x| x.to_string()).collect::<Vec<_>>().join(","),
        )?;
        Ok(())
    }

    fn write_surface(&mut self, tag: u64, lineloop: u64, isplane: bool) -> io::Result<()> {
        writeln!(self.0, "{}Surface({})={{{}}};",
            if isplane { "Plane " } else { "" }, 
            tag,
            lineloop,
        )?;
        Ok(())
    }

    fn write_physical_surface(&mut self, tag: u64, surfaces: &[u64]) -> io::Result<()> {
        writeln!(self.0, "Physical Surface(\"ps{}\")={{{}}};", tag,
            surfaces.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
        )?;
        Ok(())
    }

    fn write_surface_loop(&mut self, tag: u64, surfaces: &[i64]) -> io::Result<()> {
        writeln!(self.0, "Surface Loop({})={{{}}};",
            tag, 
            surfaces.iter().map(|&x| x.to_string()).collect::<Vec<_>>().join(","),
        )?;
        Ok(())
    }

    fn write_volume(&mut self, tag: u64, surfaceloop: u64) -> io::Result<()> {
        writeln!(self.0, "Volume({})={{{}}};", tag, surfaceloop)?;
        Ok(())
    }

    pub fn write_elem(&mut self, tag: u64, elem: &GeoElem) -> io::Result<()> {
        match elem {
            GeoElem::Point{ pos }                 => self.write_point(tag, pos)?,
            GeoElem::Line{ points }               => self.write_line(tag, points)?,
            GeoElem::LineLoop{ lines }            => self.write_line_loop(tag, lines)?,
            GeoElem::Surface{ lineloop, isplane } => self.write_surface(tag, *lineloop, *isplane)?,
            GeoElem::PhysicalSurface{ surfaces }  => self.write_physical_surface(tag, surfaces)?,
            GeoElem::SurfaceLoop{ surfaces }      => self.write_surface_loop(tag, surfaces)?,
            GeoElem::Volume{ surfaceloop }        => self.write_volume(tag, *surfaceloop)?,
        }
        Ok(())
    }

    pub fn write_geometry(&mut self, geom: &Geometry) -> io::Result<()> {
        for (&tag, elem) in 
            geom.iter(GeoElemKind::Point)
                .chain(geom.iter(GeoElemKind::Line))
                .chain(geom.iter(GeoElemKind::LineLoop))
                .chain(geom.iter(GeoElemKind::Surface))
                .chain(geom.iter(GeoElemKind::PhysicalSurface))
                .chain(geom.iter(GeoElemKind::SurfaceLoop))
                .chain(geom.iter(GeoElemKind::Volume)) {
            self.write_elem(tag, elem)?;
        }
        Ok(())
    }
}

type GeoIdxMap = indexmap::IndexMap<u64, ByAddr<Box<GeoElem>>, FnvBuildHasher>;

#[derive(Debug, Clone)]
pub struct Geometry {
    points: GeoIdxMap,
    lines: GeoIdxMap,
    lineloops: GeoIdxMap,
    surfaces: GeoIdxMap,
    physicalsurfaces: GeoIdxMap,
    surfaceloops: GeoIdxMap,
    volumes: GeoIdxMap,
}

impl Geometry {
    pub fn new() -> Self {
        Geometry{
            points: Default::default(),
            lines: Default::default(),
            lineloops: Default::default(),
            surfaces: Default::default(),
            physicalsurfaces: Default::default(),
            surfaceloops: Default::default(),
            volumes: Default::default(),
        }
    }

    pub fn is_empty(&self, kind: GeoElemKind) -> bool {
        self.apply_to_hmap(kind, |hm: &GeoIdxMap| hm.is_empty())
    }

    pub fn clear(&mut self, kind: GeoElemKind) {
        self.apply_to_hmap_mut(kind, |hm: &mut GeoIdxMap| hm.clear())
    }

    pub fn tags(&self, kind: GeoElemKind) -> impl Iterator<Item=&u64> {
        self.apply_to_hmap(kind, |hm: &GeoIdxMap| hm.keys())
    }

    pub fn elems(&self, kind: GeoElemKind) -> impl Iterator<Item=&GeoElem> {
        self.apply_to_hmap(kind, |hm: &GeoIdxMap| hm.values()).map(|x| &***x)
    }

    pub fn elems_mut(&mut self, kind: GeoElemKind) -> impl Iterator<Item=&mut GeoElem> {
        self.apply_to_hmap_mut(kind, |hm: &mut GeoIdxMap| hm.values_mut()).map(|x| &mut ***x)
    }

    pub fn iter(&self, kind: GeoElemKind) -> impl Iterator<Item=(&u64, &GeoElem)> {
        self.tags(kind).zip(self.elems(kind))
    }

    pub fn iter_mut(&mut self, kind: GeoElemKind) -> impl Iterator<Item=(&u64, &mut GeoElem)> {
        self.apply_to_hmap_mut(kind, |hm: &mut GeoIdxMap| hm.iter_mut()).map(|(k, v)| (k, &mut ***v))
    }

    pub fn get(&self, tag: u64, kind: GeoElemKind) -> Option<&GeoElem> {
        assert_ne!(tag, 0);
        self.apply_to_hmap(kind, |hm: &GeoIdxMap| hm.get(&tag).map(|x| &*x.0))
    }

    pub fn get_mut(&mut self, tag: u64, kind: GeoElemKind) -> Option<&mut GeoElem> {
        assert_ne!(tag, 0);
        self.apply_to_hmap_mut(kind, |hm: &mut GeoIdxMap| hm.get_mut(&tag).map(|x| &mut *x.0))
    }

    pub fn insert(&mut self, tag: u64, elem: Box<GeoElem>) -> Option<Box<GeoElem>> {
        assert_ne!(tag, 0);
        self.apply_to_hmap_mut(elem.kind(), |hm: &mut GeoIdxMap| hm.insert(tag, ByAddr(elem)).map(|x| x.0))
    }

    pub fn insert_from_expr(&mut self, expr: &str) -> Option<Box<GeoElem>> {
        let (tag, elem) = GeoElem::from_expr(expr)?;
        self.insert(tag, elem.into())
    }

    pub fn remove(&mut self, tag: u64, kind: GeoElemKind) -> Option<Box<GeoElem>> {
        assert_ne!(tag, 0);
        self.apply_to_hmap_mut(kind, |hm: &mut GeoIdxMap| hm.remove(&tag).map(|x| x.0))
    }

    fn apply_to_hmap<'a, R>(&'a self, kind: GeoElemKind, f: impl FnOnce(&'a GeoIdxMap) -> R) -> R {
        match kind {
            GeoElemKind::Point           => f(&self.points),
            GeoElemKind::Line            => f(&self.lines),
            GeoElemKind::LineLoop        => f(&self.lineloops),
            GeoElemKind::Surface         => f(&self.surfaces),
            GeoElemKind::PhysicalSurface => f(&self.physicalsurfaces),
            GeoElemKind::SurfaceLoop     => f(&self.surfaceloops),
            GeoElemKind::Volume          => f(&self.volumes),
        }
    }

    fn apply_to_hmap_mut<'a, R>(&'a mut self, kind: GeoElemKind, f: impl FnOnce(&'a mut GeoIdxMap) -> R) -> R {
        match kind {
            GeoElemKind::Point           => f(&mut self.points),
            GeoElemKind::Line            => f(&mut self.lines),
            GeoElemKind::LineLoop        => f(&mut self.lineloops),
            GeoElemKind::Surface         => f(&mut self.surfaces),
            GeoElemKind::PhysicalSurface => f(&mut self.physicalsurfaces),
            GeoElemKind::SurfaceLoop     => f(&mut self.surfaceloops),
            GeoElemKind::Volume          => f(&mut self.volumes),
        }
    }

    fn plane_3p_normal(&self, ptags: [u64; 3]) -> Option<na::Unit<na::Vector3<f64>>> {
        assert_ne!(ptags, [0; 3]);

        let mut poses = [na::Point3::origin(); 3];
        for (i, p) in poses.iter_mut().enumerate() {
            *p = if let GeoElem::Point{ pos } = *self.get(ptags[i], GeoElemKind::Point)? {
                pos
            } else { 
                return None; 
            };
        }

        let dirs = [poses[2] - poses[1], poses[1] - poses[0]];
        Some(na::Unit::new_normalize(dirs[0].cross(&dirs[1])))
    }

    fn plane_2l_normal(&self, ltags: [i64; 2]) -> Option<na::Unit<na::Vector3<f64>>> {
        assert_ne!(ltags, [0; 2]);
        assert_ne!(ltags, [std::i64::MIN; 2]);

        let mut ptags = [0; 4];
        for i in 0..2 {
            if let GeoElem::Line{ points } = *self.get(ltags[i].abs() as u64, GeoElemKind::Line)? {
                if ltags[i] > 0 {
                    ptags[i * 2] = points[0];
                    ptags[i * 2 + 1] = points[1];
                } else {
                    ptags[i * 2] = points[1];
                    ptags[i * 2 + 1] = points[0];
                }
                
            } else { 
                return None; 
            };
        }
        let ptags = ptags;

        let pt_eq = |x, y| ptags[x] == ptags[y];
        let collect_pt = |x, y, z| [ptags[x], ptags[y], ptags[z]];
        let ptags = 
        if pt_eq(1, 2) || pt_eq(0, 2) {
            collect_pt(0, 1, 3)
        } else if pt_eq(0, 1) || pt_eq(1, 3) {
            collect_pt(0, 2, 3)
        } else if pt_eq(0, 3) || pt_eq(2, 3) {
            collect_pt(0, 1, 2)
        } else {
            return None;
        };

        self.plane_3p_normal(ptags)
    }

    fn actual_surface_flatness_if_differs(&self, stag: u64) -> Option<FltDiff> {
        assert_ne!(stag, 0);

        let (lltag, cur_flatness) = if let GeoElem::Surface{ lineloop, isplane } = *self.get(stag, GeoElemKind::Surface)? {
            (lineloop, isplane)
        } else { 
            return None; 
        };

        let ltags = if let GeoElem::LineLoop{ ref lines } = *self.get(lltag, GeoElemKind::LineLoop)? {
            lines.clone()
        } else { 
            return None; 
        };

        let normals: Vec<_> = ltags.windows(2).map(|ls| self.plane_2l_normal([ls[0], ls[1]]).unwrap()).collect();
        let act_flatness = normals.windows(2).all(|ns| 1.0 - ns[0].dot(&ns[1]) <= f64::EPSILON);
        if cur_flatness != act_flatness {
            Some(FltDiff::Different(act_flatness))
        } else {
            Some(FltDiff::Same)
        }
    }

    pub fn correct_surface_flatness(&mut self, stag: u64) -> Option<()> {
        assert_ne!(stag, 0);

        let flatness = if let FltDiff::Different(flt) = self.actual_surface_flatness_if_differs(stag)? {
            flt
        } else {
            return Some(());
        };
        if let GeoElem::Surface{ ref mut isplane, .. } = *self.get_mut(stag, GeoElemKind::Surface).unwrap() {
            *isplane = flatness;
            Some(())
        } else {
            None
        }
    } 
}

impl From<GeoFile> for Geometry {
    fn from(mut file: GeoFile) -> Self {
        let mut geom = Self::new();
        for expr in file.read_exprs() {
            geom.insert_from_expr(&expr);
        }

        geom
    }
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
enum FltDiff {
    Same,
    Different(bool),
}

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum GeoElemKind {
    Point,
    Line,
    LineLoop,
    Surface,
    PhysicalSurface,
    SurfaceLoop,
    Volume,
}

#[derive(Debug, Clone)]
pub enum GeoElem {
    Point{
        pos: na::Point3<f64>,
    },
    Line{
        points: [u64; 2],
    },
    LineLoop{
        lines: Vec<i64>,
    },
    Surface{
        lineloop: u64,
        isplane: bool,
    },
    PhysicalSurface{
        surfaces: Vec<u64>,
    },
    SurfaceLoop{
        surfaces: Vec<i64>,
    },
    Volume{
        surfaceloop: u64,
    },
}

impl GeoElem {
    pub fn kind(&self) -> GeoElemKind {
        match self {
            GeoElem::Point{ .. }           => GeoElemKind::Point,
            GeoElem::Line{ .. }            => GeoElemKind::Line,
            GeoElem::LineLoop{ .. }        => GeoElemKind::LineLoop,
            GeoElem::Surface{ .. }         => GeoElemKind::Surface,
            GeoElem::PhysicalSurface{ .. } => GeoElemKind::PhysicalSurface,
            GeoElem::SurfaceLoop{ .. }     => GeoElemKind::SurfaceLoop,
            GeoElem::Volume{ .. }          => GeoElemKind::Volume,
        }
    } 

    pub fn from_expr(expr: &str) -> Option<(u64, Self)> {
        let tag: u64 = str_between(expr, "(", ")")?.parse().ok()?;
        let other: Vec<&str> = str_between(expr, "{", "}")?.split(",").map(|x| x.trim()).collect();

        if expr.starts_with("Point") {
            let parse_res: Vec<_> = other.into_iter().map(|x| x.parse::<f64>()).collect::<Result<_, _>>().ok()?;
            let pos = na::Point3::from_slice(&parse_res);
            Some((tag, GeoElem::Point{ pos }))

        } else if expr.starts_with("Line Loop") {
            let lines: Vec<_> = other.into_iter().map(|x| x.parse::<i64>()).collect::<Result<_, _>>().ok()?;
            Some((tag, GeoElem::LineLoop{ lines }))

        } else if expr.starts_with("Line") {
            let points: Vec<_> = other.into_iter().map(|x| x.parse::<u64>()).collect::<Result<_, _>>().ok()?;
            Some((tag, GeoElem::Line{ points: [points[0], points[1]] }))

        } else if expr.starts_with("Plane Surface") {
            if other.len() != 1 {
                return None;
            }
            let lineloop = other[0].parse::<u64>().ok()?;
            let isplane = true;
            Some((tag, GeoElem::Surface{ lineloop, isplane }))

        } else if expr.starts_with("Physical Surface") {
            let surfaces: Vec<_> = other.into_iter().map(|x| x.parse::<u64>()).collect::<Result<_, _>>().ok()?;
            Some((tag, GeoElem::PhysicalSurface{ surfaces }))

        } else if expr.starts_with("Surface Loop") {
            let surfaces: Vec<_> = other.into_iter().map(|x| x.parse::<i64>()).collect::<Result<_, _>>().ok()?;
            Some((tag, GeoElem::SurfaceLoop{ surfaces }))

        } else if expr.starts_with("Surface") {
            if other.len() != 1 {
                return None;
            }
            let lineloop = other[0].parse::<u64>().ok()?;
            let isplane = false;
            Some((tag, GeoElem::Surface{ lineloop, isplane }))

        } else if expr.starts_with("Volume") {
            if other.len() != 1 {
                return None;
            }
            let surfaceloop = other[0].parse::<u64>().ok()?;
            Some((tag, GeoElem::Volume{ surfaceloop }))

        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let file = GeoFile::open("test.geo").unwrap();
        let mut geom = Geometry::from(file);
        geom.clear(GeoElemKind::PhysicalSurface);
        let stags: Vec<u64> = geom.tags(GeoElemKind::Surface).map(|x| *x).collect();
        for stag in stags {
            geom.correct_surface_flatness(stag).unwrap();
        }
        let mut file = OccFile::create("test-occ.geo").unwrap();
        file.write_geometry(&geom).unwrap();
    }
}
