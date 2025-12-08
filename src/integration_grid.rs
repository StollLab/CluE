use crate::space_3d::Vector3D;
use crate::clue_errors::*;
use crate::io::FromTOMLString;
use crate::io::write_data;
use crate::physical_constants::PI;

use lebedev_laikov;
use rand_chacha::ChaCha20Rng;

use serde::{Serialize,Deserialize};
use std::fmt;

/// `IntegrationGrid` is a struct for use in 
/// [quadrature](https://en.wikipedia.org/wiki/Quadrature_(mathematics)).
#[derive(Debug, Clone,PartialEq,Serialize,Deserialize)]
pub struct IntegrationGrid{
  pub dim: usize,
  pub points: Vec::<f64>,
  pub weights: Vec::<f64>,
}

impl FromTOMLString for IntegrationGrid{
  fn from_toml_string(toml_str: &str) -> Result<Self,CluEError>{
    let decoded: Result<IntegrationGrid,_> = toml::from_str(toml_str);
    match decoded {
      Ok(out) => Ok(out),
      Err(err) => Err(CluEError::CannotReadTOML( format!("{}",err) )),
    }
  }
}
impl fmt::Display for IntegrationGrid{
   fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
     let toml = toml::to_string(self).unwrap();
     write!(f,"{}",toml)
   }
} 

impl IntegrationGrid{

  /// This function generates a 
  /// [Lebedev](https://en.wikipedia.org/wiki/Lebedev_quadrature)
  /// quadrature grid as an `IntegrationGrid`.
  /// It uses the
  /// [lebedev_laikov crate](https://github.com/Rufflewind/lebedev_laikov)
  /// by [Rufflewind](https://github.com/Rufflewind).
  /// # Example
  /// ```
  ///# use clue_oxide::integration_grid::IntegrationGrid; 
  /// let n = 26;
  /// let grid = IntegrationGrid::lebedev(n);
  /// ```
  ///
  /// # Panic
  /// `lebedev()` will paninc if
  /// `n` ∉  \{6,   14,   26,   38,   50,   74,   86,  110,  146,  170,
  ///        194,  230,  266,  302,  350,  434,  590,  770,  974, 1202,
  ///       1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802,
  ///       5294, 5810\}.
  ///
  pub fn lebedev(n: usize) -> Result<Self,CluEError>{

    const LEBEDEVGRIDSIZES: [usize;32] = [
               6,   14,   26,   38,   50,   74,   86,  110,  146,  170, 
             194,  230,  266,  302,  350,  434,  590,  770,  974, 1202, 
            1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 
            5294, 5810];

    if !LEBEDEVGRIDSIZES.contains(&n){
      return Err(CluEError::NotALebedevGrid(n));
    }

    let (x, y, z, weights) = lebedev_laikov::ld_vecs(n);

    let dim = 3;
    let n_elements = n*dim;

    let mut points = Vec::<f64>::with_capacity(n_elements);

    for ((xi,yi),zi) in x.iter().zip(y.iter()).zip(z.iter()){
      points.push(*xi);
      points.push(*yi);
      points.push(*zi);

    }

    Ok(IntegrationGrid{dim,points,weights})

  }

  //----------------------------------------------------------------------------
  /// This function generates a single point grid one unit in the z-direction.
  pub fn z_3d() -> Self{
    let dim = 3;
    let points = vec![0.0, 0.0, 1.0];
    let weights = vec![1.0];
    IntegrationGrid{dim,points,weights}
  }
  //----------------------------------------------------------------------------
  /// This function generates a new `dim`-dimensional grid with no points.
  pub fn new(dim: usize) -> Self {
    let points = Vec::<f64>::new();
    let weights = Vec::<f64>::new();

    IntegrationGrid{dim,points,weights}
  }
  //----------------------------------------------------------------------------
  /// This function generations a grid with `num_point` points randomly 
  /// distributed over the unit sphere.
  pub fn random_unit_sphere(num_points: usize,rng: &mut ChaCha20Rng) -> Self
  {
    let dim = 3;
    let mut points = Vec::<f64>::with_capacity(dim*num_points);
    let w = 1.0/(num_points as f64);
    let weights = (0..num_points).map(|_| w).collect::<Vec::<f64>>();

    for _ii in 0..num_points{
      let r = Vector3D::random_direction(rng);
      points.push(r.x());
      points.push(r.y());
      points.push(r.z());
    }

    IntegrationGrid{dim,points,weights}
  }
  //----------------------------------------------------------------------------
  /// This function scales all the points by a `scale_factor`.
  pub fn scale(&self, scale_factor: f64) -> Self{

    let points = self.points.iter().map(|x| x*scale_factor)
      .collect::<Vec::<f64>>();
    IntegrationGrid{dim: self.dim ,points,weights: self.weights.clone()}

  }
  //----------------------------------------------------------------------------
  /// This function pushes a new point to a grid.
  pub fn push(&mut self, points: Vec::<f64>,weight: f64)
    -> Result<(),CluEError>
  {
    if points.len() != self.dim{
      return Err(CluEError::CannotAddPointToGrid(points.len(), self.dim));
    }

    for x in points.into_iter(){
      self.points.push(x);
    }

    self.weights.push(weight);

    Ok(())
  }
  //----------------------------------------------------------------------------
  /// This function returns the number of points in a grid.
  pub fn len(&self) -> usize { self.weights.len() }
  //----------------------------------------------------------------------------
  /// This function returns the dimensionality of the space that the grid
  /// is imbeded in. 
  pub fn dim(&self) -> usize { self.dim }
  //----------------------------------------------------------------------------
  /// This function checks if the grid is empty.
  pub fn is_empty(&self) -> bool { self.weights.is_empty() }
  //----------------------------------------------------------------------------
  /// This function returns the weight of point `index`.
  pub fn weight(&self, index: usize) -> f64{ 
    self.weights[index]
  }
  //----------------------------------------------------------------------------
  /// This function returns the x value of point `index`.
  pub fn x(&self, index: usize) -> f64{ 
    let idx = self.dim*index;
    self.points[idx]
  }
  //----------------------------------------------------------------------------
  /// This function returns the y value of point `index`.
  /// The function will panic if the space is not at least 2D.
  pub fn y(&self, index: usize) -> f64{ 
    assert!(self.dim >= 2);
    let idx = self.dim*index;
    self.points[idx+1]
  }
  //----------------------------------------------------------------------------
  /// This function returns the z value of point `index`.
  /// The function will panic if the space is not at least 3D.
  pub fn z(&self, index: usize) -> f64{ 
    assert!(self.dim >= 3);
    let idx = self.dim*index;
    self.points[idx+2]
  }
  //----------------------------------------------------------------------------
  /// This function returns the coordinates of point `index`.
  /// The function will err if the space is not 3D.
  pub fn xyz(&self, index: usize) -> Result<Vector3D,CluEError>{
    if self.dim != 3{
      return Err(CluEError::NotA3DVector(self.dim));
    }
    Ok(Vector3D::from( [ self.x(index), self.y(index), self.z(index) ] ))
  }
  //----------------------------------------------------------------------------
  /// This function returns the centroid of the grid.
  pub fn mean(&self) -> Vec::<f64>{
    let mut r: Vec::<f64> = (0..self.dim).map(|_| 0.0).collect();

    for (ii,x) in self.points.iter().enumerate(){
      let idx = ii%self.dim;
      r[idx] += x;
    }

    let normalizing_factor = (self.dim as f64)/(self.points.len() as f64);
    r = r.iter().map(|x| normalizing_factor*x).collect();

    r
  }
  //----------------------------------------------------------------------------
  /// This function translates all the points in the grid by `vector`.
  /// The function will panic if the vector is of a different dimension than
  /// the grid.
  pub fn translate_mut(&mut self, vector: &[f64]){
    assert_eq!(vector.len(), self.dim);

    for (ii,x) in self.points.iter_mut().enumerate(){
      let idx = ii%self.dim;
      *x += vector[idx];
    }

  }
  //----------------------------------------------------------------------------
  /// This function reads a grid from a csv file. 
  pub fn read_from_csv(filename: &str) -> Result<Self,CluEError>{

    let Ok(mut rdr) = csv::Reader::from_path(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()));
    };

    let mut num_rows = 0;
    let mut num_cols = 0;
    
    for result in rdr.records() {
      if let Ok(str_rec) = result{
        num_cols = str_rec.len();
      }else{
        return Err(CluEError::CannotOpenFile(filename.to_string()));
      };
      num_rows += 1;
    } 

    if num_cols <= 1{
      return Err(CluEError::CannotReadGrid(filename.to_string()));
    }
    let dim = num_cols - 1;

    let mut points = Vec::<f64>::with_capacity(dim*num_rows);
    let mut weights = Vec::<f64>::with_capacity(num_rows);
    
    let Ok(mut rdr) = csv::Reader::from_path(filename) else{
      return Err(CluEError::CannotOpenFile(filename.to_string()));
    };

    for result in rdr.records() {
      let Ok(record) = result else{
        return Err(CluEError::CannotOpenFile(filename.to_string()));
      };

      for (icol,entry) in record.iter().enumerate(){
        let Ok(v) = entry.parse::<f64>() else{
          return Err(CluEError::CannotOpenFile(filename.to_string()));
        };
        if icol < dim{
          points.push(v);
        }else{
          weights.push(v);
        }
      }
    }

    if points.len() != dim*num_rows || weights.len() != num_rows{
      return Err(CluEError::CannotReadGrid(filename.to_string()));
    }

    Ok(IntegrationGrid{dim,points,weights})

  }
  //----------------------------------------------------------------------------
  pub fn to_csv(&self,filename: &str) -> Result<(),CluEError>
  {
    let mut  headers = Vec::<String>::with_capacity(self.dim+1);
    for idx in 1..=self.dim{
      headers.push(format!("x{}",idx)); 
    }
    headers.push("w".to_string()); 

    let mut data = Vec::<Vec::<f64>>::with_capacity(self.dim+1);
    for _ in 0..self.dim+1{
      data.push(Vec::<f64>::with_capacity(self.len()));
    }

    for ii in  0..self.len(){
      for ix in 0..self.dim{
        data[ix].push(self.points[ii + ix]);
      }
      data[self.dim].push(self.weights[ii]);
    }

    write_data(&data,filename,headers)

  }
  //----------------------------------------------------------------------------
  pub fn to_polar_csv(&self,filename: &str) -> Result<(),CluEError>
  {
    assert_eq!(self.dim,3);

    let headers = vec![
      "r".to_string(),
      "theta".to_string(),
      "phi".to_string(),
      "weight".to_string(),
      "theta_degrees".to_string(),
      "phi_degrees".to_string(),
    ];
      

    let mut data = Vec::<Vec::<f64>>::with_capacity(5);
    for _ in 0..6{
      data.push(Vec::<f64>::with_capacity(self.len()));
    }

    for ii in  0..self.len(){
      let x = self.x(ii);
      let y = self.y(ii);
      let z = self.z(ii);
      let r = (x*x + y*y + z*z).sqrt();

      let theta = (z/r).acos();
      let rho = (x*x + y*y).sqrt();
      let mut phi = (x/rho).acos();
      if y < 0.0{
        phi = -phi;
      }

      data[0].push(r);
      data[1].push(theta);
      data[2].push(phi);
      data[3].push(self.weights[ii]);
      data[4].push(theta*180.0/PI);
      data[5].push(phi*180.0/PI);
    }

    write_data(&data,filename,headers)

  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  /// This function removes the negative hemisphere of a grid.
  pub fn remove_3d_hemisphere(&self) -> IntegrationGrid{
    assert_eq!(self.dim,3);

    let mut points = Vec::<f64>::with_capacity(self.points.len()/2);
    let mut weights = Vec::<f64>::with_capacity(self.len()/2);

    let thr = 1e-12;
    for ii in 0..self.len(){
      let x = self.x(ii);
      let y = self.y(ii);
      let z = self.z(ii); 
      
      if (z < 0.0) || (z < thr && y < 0.0) || (z < thr && y < thr && x < 0.0){ 
        continue; 
      }

      points.push(x);
      points.push(y);
      points.push(z);
      weights.push(self.weight(ii));
    }

    let norm: f64 = weights.iter().sum();

    IntegrationGrid{
      dim: self.dim,
      points,
      weights: weights.iter().map(|x| x/norm).collect()
    }
  
  }
  //----------------------------------------------------------------------------
 
}



#[cfg(test)]
mod tests{
  use super::*;

  //----------------------------------------------------------------------------
  #[test]
  fn test_read_from_csv(){
    let grid = IntegrationGrid::read_from_csv("assets/xyzw.csv").unwrap();
    assert_eq!(grid.dim, 3);
    assert_eq!(grid.points, vec![1.0,0.0,0.0,
                                 0.0,1.0,0.0,
                                 0.0,0.0,1.0]);
    assert_eq!(grid.weights, vec![0.3333,0.3333,0.3333]);

  }
  //----------------------------------------------------------------------------
  #[test]
  fn test_lebedev(){
    let grid = IntegrationGrid::lebedev(6).unwrap();
    let norm: f64 = grid.weights.iter().sum();
    assert!( (norm- 1.0).abs() <1e12);

    assert_eq!(grid.points, vec![
        1.0, 0.0, 0.0,
        -1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, 0.0, 1.0,
        0.0, 0.0, -1.0,
    ]);

    assert_eq!(grid.mean(), vec![0.0,0.0,0.0]);
  } 
  //----------------------------------------------------------------------------
  #[test]
  fn test_remove_3d_hemisphere(){
    let grid = IntegrationGrid::lebedev(6).expect("not a Lebedev grid")
      .remove_3d_hemisphere();
    let norm: f64 = grid.weights.iter().sum();
    assert!( (norm- 1.0).abs() <1e12);
    assert_eq!(grid.points, vec![
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    ]);

    assert_eq!(grid.mean(), vec![1.0/3.0,1.0/3.0,1.0/3.0]);

  } 
  //----------------------------------------------------------------------------
  #[test]
  fn test_translate_mut(){
    let mut grid = IntegrationGrid::new(3);
    let w = 1.0/6.0;
    grid.push(vec![1.0, 20.0, 200.0],w).unwrap();
    grid.push(vec![3.0, 20.0, 200.0],w).unwrap();
    grid.push(vec![2.0, 10.0, 200.0],w).unwrap();
    grid.push(vec![2.0, 30.0, 200.0],w).unwrap();
    grid.push(vec![2.0, 20.0, 100.0],w).unwrap();
    grid.push(vec![2.0, 20.0, 300.0],w).unwrap();
    assert_eq!(grid.mean(), vec![2.0, 20.0, 200.0]);
    let r = vec![-2.0,-20.0,-200.0];
    grid.translate_mut(&r);
    assert_eq!(grid.mean(), vec![0.0, 0.0, 0.0]);
  }
}






