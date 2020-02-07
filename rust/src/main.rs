
use plotters::prelude::*;
use rand::Rng;


macro_rules! rng {
    () => {
    {
        let mut rng = rand::thread_rng();
        let a = 2.0*rng.gen::<f32>() - 1.0;
        a
    }
    };
}

const G: f32 = 1.0;


trait Calculate {
    fn next_step(&self, h: f32, t: f32, par: &Vec<Particle>, sys: System) -> Vec<Particle>;
}

trait EquationsOfMotion {
    fn dq(&self, i: usize, t: f32, q: &Vec<Vector3>, p: &Vec<Vector3>, m: &Vec<f32>) -> Vector3;
    fn dp(&self, i: usize, t: f32, q: &Vec<Vector3>, p: &Vec<Vector3>, m: &Vec<f32>) -> Vector3;
}

enum Method {
    Leapfrog,
    RungeKutta,
    Euler
}

#[derive(Copy, Clone)]
enum System {
    Gravitational,
}

impl Calculate for Method {
    fn next_step(&self, h: f32, t: f32, par: &Vec<Particle>, sys: System) -> Vec<Particle> {
        match self {
            Method::Leapfrog => {
                let mut p0: Vec<Vector3> = Vec::new();
                let mut q0: Vec<Vector3> = Vec::new();
                let mut m: Vec<f32> = Vec::new();
                for part in par.iter() {
                    p0.push(part.p);
                    q0.push(part.q);
                    m.push(part.m);
                } 

                let mut p12: Vec<Vector3> = Vec::new();
                for i in 0..par.len() {
                    p12.push(Vector3::addition(p0[i], sys.dp(i, t, &q0, &p0, &m).scale(h*0.5)));
                }

                let mut q1: Vec<Vector3> = Vec::new();
                for i in 0..par.len() {
                    q1.push(Vector3::addition(q0[i], sys.dq(i, t, &q0, &p12, &m).scale(h)));
                }

                let mut p1: Vec<Vector3> = Vec::new();
                for i in 0..par.len() {
                    p1.push(Vector3::addition(p12[i], sys.dp(i, t, &q1, &p12, &m).scale(h*0.5)));
                }

                let mut par1: Vec<Particle> = Vec::new();
                for i in 0..par.len() {
                    par1.push(Particle::new(q1[i], p1[i], m[i]));
                }

                return par1;
            },
            _ => {
                return Vec::new();
            }
        }
    }
}

impl EquationsOfMotion for System {
    fn dq(&self, i: usize, t: f32, q: &Vec<Vector3>, p: &Vec<Vector3>, m: &Vec<f32>) -> Vector3 {
        match self {
            System::Gravitational => {
                return p[i];
            },
            _ => {
                return p[i];
            }
        }
    }

    fn dp(&self, i: usize, t: f32, q: &Vec<Vector3>, p: &Vec<Vector3>, m: &Vec<f32>) -> Vector3 {
        match self {
            System::Gravitational => {
                let mut a = Vector3::new(0.0, 0.0, 0.0);
                let qi = q[i];
                for j in 0..q.len() {
                    if j != i {
                        let r = Vector3::substraction(q[j], qi);
                        let l = r.lenght();
                        let c = G*m[j]/(l*l);
                        let aj = r.unit_vector().scale(c);
                        a = Vector3::addition(a, aj);
                    }
                }
                return a;
            },
            _ => {
                return Vector3::new(0.0, 0.0, 0.0);
            }
        }
    }
}


#[derive(Copy, Clone)]
struct Vector3 {
    x: f32,
    y: f32,
    z: f32
}

impl Vector3 {
    fn new(x: f32, y: f32, z: f32) -> Vector3 {
        return Vector3{x, y, z};
    }

    fn lenght(&self) -> f32 {
        let l = self.x * self.x + self.y * self.y + self.z * self.z;
        return l.sqrt();
    }

    fn unit_vector(&self) -> Vector3 {
        let l = self.lenght();
        let u = Vector3::new(self.x/l, self.y/l, self.z/l);
        return u;
    }

    fn addition(a: Vector3, b: Vector3) -> Vector3 {
        let v = Vector3::new(a.x+b.x, a.y+b.y, a.z+b.z);
        return v;
    }
    
    fn substraction(a: Vector3, b: Vector3) -> Vector3 {
        let v = Vector3::new(a.x-b.x, a.y-b.y, a.z-b.z);
        return v;
    }

    fn scale(&self, c: f32) -> Vector3 {
        return Vector3::new(c*self.x, c*self.y, c*self.z);
    }
}

struct Particle {
    q: Vector3,
    p: Vector3,
    m: f32
}

impl Particle {
    fn new(q: Vector3, p: Vector3, m: f32) -> Particle {
        return Particle{q, p, m};
    }

    fn center_of_mass(par: &Vec<Particle>) -> Vector3 {
        let mut cm = Vector3::new(0.0, 0.0, 0.0);
        let mut m = 0.0;
        for i in 0..par.len() {
            m += par[i].m;
            cm = Vector3::addition(cm, par[i].q.scale(par[i].m));
        }
        return cm.scale(1.0/m);
    }
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::gif("test.gif", (1024, 768), 17)?.into_drawing_area();
    
    let n = 400;

    let mut par: Vec<Particle> = vec![];

    let mut par: Vec<Particle> = vec![
        Particle::new(Vector3::new(0.0, -3.0, 0.0), Vector3::new(-0.5, 0.0, 0.0), 5.0),
        Particle::new(Vector3::new(0.0, 3.0, 0.0), Vector3::new(0.5, 0.0, 0.0), 5.0),
        //Particle::new(Vector3::new(-1.0, 0.0, 0.0), Vector3::new(0.5, -0.5, 0.0), 1.0)
    ]; 

    //par.push(Particle::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0), 100.0));

    //for _ in 0..n {
    //    par.push(Particle::new(Vector3::new(rng![]*4.0, rng![]*4.0, rng![]*4.0), Vector3::new(rng![], rng![], rng![]), 1.0));
    //}

    let mut t = 0.0;
    let h = 0.005;
    let steps = 20000;

    let afa: u32 = 20;


    let meth: Method = Method::Leapfrog;
    let sys: System = System::Gravitational;

    for i in 0..steps {

        par = meth.next_step(h, t, &par, sys);    
        t = t + h;

        if i%afa == 0 {

            let cm: Vector3 = Particle::center_of_mass(&par);
            // let cm: Vector3 = par[0].q;

            let title = format!("step: {}", i);

            root.fill(&WHITE)?;
            let mut chart = ChartBuilder::on(&root)
                .x_label_area_size(40)
                .y_label_area_size(40)
                .margin(30)
                .caption(title, ("sans-serif", 25).into_font())
                .build_ranged((cm.x-5.0)..(cm.x+5.0), (cm.y-5.0)..(cm.y+5.0))?;
        
            chart.configure_mesh().line_style_2(&WHITE).draw()?;
            
            let mut points: Vec<(f32, f32)> = Vec::new();
            for pi in par.iter() {
                points.push((pi.q.x, pi.q.y));
            }

            chart.draw_series(
                points
                    .iter()
                    .map(|(x, y)| Circle::new((*x, *y), 4, BLUE.filled())),
            )?;
            root.present()?;
        }

        if i%(steps/100) == 0 {
            println!("complete: {}%", (i*100)/steps);
        }

    }
    Ok(())
}

