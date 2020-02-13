
use plotters::prelude::*;
use rand::{SeedableRng, Rng, rngs::StdRng};
use std::time::{Duration, SystemTime};


macro_rules! rng {
    ($rng:expr) => {
    {   
        let a = 2.0*$rng.gen::<f32>() - 1.0;
        a
    }
    };
}

const G: f32 = 1.0;
const COL: f32 = 0.04;

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

impl Method {
    fn collision_check(par: &Vec<Particle>, col: f32, del: f32, cm: Vector3) -> (Vec<Particle>, Vec<u32>) {
        let mut par1: Vec<Particle> = vec![];
        let mut destroyed: Vec<u32> = vec![];
        let mut done: Vec<usize> = vec![];
        let mut ind = 0;
        'check: for i in 0..par.len() {
            if Vector3::substraction(cm, par[i].q).lenght() > del {
                destroyed.push(par[i].id);
                continue 'check;
            }
            for u in done.iter() {
                if i == *u {
                    destroyed.push(par[i].id);
                    continue 'check;
                }
            }
            for j in i+1..par.len() {
                let r = Vector3::substraction(par[i].q, par[j].q);
                if r.lenght() < col {
                    let m = par[i].m + par[j].m;
                    let q = Vector3::addition(par[i].q.scale(par[i].m/m), par[j].q.scale(par[j].m/m));
                    let p = Vector3::addition(par[i].p.scale(par[i].m/m), par[j].p.scale(par[j].m/m));
                    let np = Particle::new(q, p, m, ind);
                    par1.push(np);
                    ind += 1;
                    done.push(j);
                    continue 'check;
                }
            }
            par1.push(Particle::new(par[i].q, par[i].p, par[i].m, ind));
            ind += 1;
        }

        return (par1, destroyed);
    }
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
                    par1.push(Particle::new(q1[i], p1[i], m[i], par[i].id));
                }

                return par1;
            },
            Method::RungeKutta => {
                let mut p0: Vec<Vector3> = Vec::new();
                let mut q0: Vec<Vector3> = Vec::new();
                let mut m: Vec<f32> = Vec::new();
                for part in par.iter() {
                    p0.push(part.p);
                    q0.push(part.q);
                    m.push(part.m);
                }

                let mut k1: Vec<Vector3> = Vec::new();
                let mut l1: Vec<Vector3> = Vec::new();
                for i in 0..par.len() {
                    k1.push(sys.dq(i, t, &q0, &p0, &m).scale(h));
                    l1.push(sys.dp(i, t, &q0, &p0, &m).scale(h));
                }

                let mut k2: Vec<Vector3> = Vec::new();
                let mut l2: Vec<Vector3> = Vec::new();
                let mut q01: Vec<Vector3> = Vec::new();
                let mut p01: Vec<Vector3> = Vec::new();
                for i in 0..par.len() {
                    q01.push(Vector3::addition(q0[i], k1[i].scale(0.5)));
                    p01.push(Vector3::addition(p0[i], l1[i].scale(0.5)));
                }
                for i in 0..par.len() {
                    k2.push(sys.dq(i, t+0.5*h, &q01, &p01, &m).scale(h));
                    l2.push(sys.dp(i, t+0.5*h, &q01, &p01, &m).scale(h));
                }

                let mut k3: Vec<Vector3> = Vec::new();
                let mut l3: Vec<Vector3> = Vec::new();
                let mut q02: Vec<Vector3> = Vec::new();
                let mut p02: Vec<Vector3> = Vec::new();
                for i in 0..par.len() {
                    q02.push(Vector3::addition(q0[i], k2[i].scale(0.5)));
                    p02.push(Vector3::addition(p0[i], l2[i].scale(0.5)));
                }
                for i in 0..par.len() {
                    k3.push(sys.dq(i, t+0.5*h, &q02, &p02, &m).scale(h));
                    l3.push(sys.dp(i, t+0.5*h, &q02, &p02, &m).scale(h));
                }
                
                let mut k4: Vec<Vector3> = Vec::new();
                let mut l4: Vec<Vector3> = Vec::new();
                let mut q03: Vec<Vector3> = Vec::new();
                let mut p03: Vec<Vector3> = Vec::new();
                for i in 0..par.len() {
                    q03.push(Vector3::addition(q0[i], k3[i]));
                    p03.push(Vector3::addition(p0[i], k3[i]));
                }
                for i in 0..par.len() {
                    k4.push(sys.dq(i, t+h, &q03, &p03, &m).scale(h));
                    l4.push(sys.dp(i, t+h, &q03, &p03, &m).scale(h));
                }


                let mut q1: Vec<Vector3> = Vec::new();
                let mut p1: Vec<Vector3> = Vec::new();
                for i in 0..par.len() {
                    let tmp_q = Vector3::addition_n(vec![k1[i], k2[i].scale(2.0), k3[i].scale(2.0), k4[i]]).scale(1.0/6.0);
                    let tmp_p = Vector3::addition_n(vec![l1[i], l2[i].scale(2.0), l3[i].scale(2.0), l4[i]]).scale(1.0/6.0);

                    q1.push(Vector3::addition(q0[i], tmp_q));
                    p1.push(Vector3::addition(p0[i], tmp_p));
                }

                let mut par1: Vec<Particle> = Vec::new();
                for i in 0..par.len() {
                    par1.push(Particle::new(q1[i], p1[i], m[i], par[i].id));
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
    
    fn addition_n(a: Vec<Vector3>) -> Vector3 { 
        let mut v = Vector3::new(0.0, 0.0, 0.0);
        for u in a.iter() {
            v = Vector3::addition(v, *u);
        }
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

#[derive(Copy, Clone)]
struct Particle {
    q: Vector3,
    p: Vector3,
    m: f32,
    id: u32
}

impl Particle {
    fn new(q: Vector3, p: Vector3, m: f32, id: u32) -> Particle {
        return Particle{q: q, p: p, m: m, id: id};
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

fn hash(n: f32) -> u32 {
    let un = (n.abs()*10000.0) as u32;
    return 7 + 17*un + 33*un*un + 23*un*un*un + 13*un*un*un*un;
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    
    // time
    let time0 = SystemTime::now();

    let mut rend_t = 0;
    let mut calc_t = 0;

    // simulation
    let meth1: Method = Method::Leapfrog;
    let meth2: Method = Method::RungeKutta;
    let sys: System = System::Gravitational;
    
    let mut n = 5000;
    let speed_mult: f32 = 2.0;
    let spawn_area: f32 = 15.0;
    let mut rng = StdRng::from_seed([5; 32]);

    let mut t = 0.0;
    let h = 0.001;
    let steps = 50000;

    let collisions: bool = true;
    let collision_dis: f32 = 0.05;
    let deletion_dis: f32 = 100.0;

    // rendering
    let line_density: u32 = 10;
    let border: f32 = 15.0;
    let afa: u32 = 50;
    let dot_size: f32 = 2.0;
    let line: bool = false;
    let line_len: usize = 50;

    let root = BitMapBackend::gif("test.gif", (1000, 1000), 50)?.into_drawing_area();


    let mut par1: Vec<Particle> = vec![];
    //let mut par2: Vec<Particle> = vec![];

    /*
    let mut par: Vec<Particle> = vec![
        Particle::new(Vector3::new(0.0, -1.3, 0.0), Vector3::new(0.5, -0.25, 0.0), 5.0, 1),
        Particle::new(Vector3::new(0.0, 1.3, 0.0), Vector3::new(-0.5, -0.25, 0.0), 5.0, 2),
        Particle::new(Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.5, 0.0), 1.0, 3),
        //Particle::new(Vector3::new(-2.0, 0.0, 0.0), Vector3::new(0.0, 0.5, 0.0), 1.0, 4)
    ]; 
    for p in par.iter() {
        par1.push(*p);
        par2.push(*p);
    }*/

    //par.push(Particle::new(Vector3::new(0.0, 0.0, 0.0), Vector3::new(0.0, 0.0, 0.0), 100.0, 100));
    
    for i in 0..n {
        let p = Particle::new(
            Vector3::new(rng![rng]*spawn_area, rng![rng]*spawn_area, rng![rng]*spawn_area), 
            Vector3::new(rng![rng]*speed_mult, rng![rng]*speed_mult, rng![rng]*speed_mult), 
            0.2,
            i
        );
        par1.push(p);
        //par2.push(p);
    }

    let mut lines1 = vec![];
    //let mut lines2 = vec![];

    if line {
        for i in 0usize..(n as usize) {
            lines1.push(vec![par1[i].q]);
            //lines2.push(vec![par2[i].q]);
        }
    }

    'main: for i in 0..steps {

        let mut step_t = SystemTime::now();
        let mut dest: Vec<u32> = vec![];
        let cm: Vector3 = Particle::center_of_mass(&par1);

        par1 = meth1.next_step(h, t, &par1, sys);
        if collisions {
            let tmp = Method::collision_check(&par1, collision_dis, deletion_dis, cm);
            par1 = tmp.0;
            dest = tmp.1;
            dest.sort();
            dest.reverse();
            if line {

                for id in dest.iter() {
                    println!("{}   {}", lines1.len(), *id);
                    lines1.remove(*id as usize);
                }
            }
        }

        n = par1.len() as u32;
        //par2 = meth2.next_step(h, t, &par2, sys);    
        t = t + h;

        calc_t += step_t.elapsed()?.as_millis();
        
        if i%line_density == 0 {
            if line {
                'line_gen: for j in 0usize..(n as usize) {
                    lines1[j].push(par1[j].q);
                    if lines1[j].len() > line_len {
                        lines1[j].remove(0);
                    }
                    //lines2[j].push(par2[j].q);
                }
            }
        }

        if i%afa == 0 {
            step_t = SystemTime::now();
            
            // let cm: Vector3 = Vector3::new(0.0, 0.0, 0.0);
            let cm: Vector3 = Particle::center_of_mass(&par1);
            // let cm: Vector3 = par1[0].q;

            let title = format!("step: {}", i);

            root.fill(&WHITE)?;
            let mut chart = ChartBuilder::on(&root)
                .x_label_area_size(0)
                .y_label_area_size(0)
                .margin(0)
                //.caption(title, ("sans-serif", 25).into_font())
                .build_ranged((cm.x-border)..(cm.x+border), (cm.y-border)..(cm.y+border))?;
        
            //chart.configure_mesh().line_style_2(&WHITE).draw()?;
            
            let mut points1: Vec<(f32, f32, f32)> = Vec::new();
            'point_loop: for pi in par1.iter() {
                if (cm.x-pi.q.x.abs()) > border || (cm.y-pi.q.y.abs()) > border {
                    continue 'point_loop;
                }
                points1.push((pi.q.x, pi.q.y, pi.m));
            }
            
            /*let mut points2: Vec<(f32, f32)> = Vec::new();
            for pi in par2.iter() {
                points2.push((pi.q.x, pi.q.y));
            }*/
            
            chart.draw_series(
                points1
                    .iter()
                    .map(|(x, y, m)| Circle::new((*x, *y), ((dot_size + *m).sqrt() as i32), BLUE.filled())),
            )?;
            //.label("Runge-Kutta")
            //.legend(|(x, y)| Circle::new((x, y), 8, BLUE.filled()));
            
            if line {
                'line_loop: for j in 0usize..(n as usize) {
                    let mut lin: Vec<(f32, f32)> = Vec::new();
                    for point in lines1[j].iter() {
                        lin.push((point.x, point.y));
                    }
                    chart.draw_series(LineSeries::new(
                        lin
                            .iter()
                            .map(|(x, y)| (*x, *y)), &BLUE
                    ))?;
                }
            }

            /*chart.draw_series(
                points2
                    .iter()
                    .map(|(x, y)| Circle::new((*x, *y), dot_size, RED.filled())),
            )?
            .label("Leapfrog")
            .legend(|(x, y)| Circle::new((x, y), 8, RED.filled()));
            
            if line {
                for j in 0..n {
                    let mut lin: Vec<(f32, f32)> = Vec::new();
                    for point in lines2[j].iter() {
                        lin.push((point.x, point.y));
                    }
                    chart.draw_series(LineSeries::new(
                        lin
                            .iter()
                            .map(|(x, y)| (*x, *y)), &RED
                    ))?;
                }
            }*/

            chart.configure_series_labels().position(SeriesLabelPosition::UpperRight).label_font(("sans-serif", 40).into_font()).draw()?;
            
            root.present()?;

            rend_t += step_t.elapsed()?.as_millis();

        }
        

        if i%(steps/500) == 0 {
            println!("complete: {}% \nparticles: {} \ncalculating: {:.2} s\nrendering:   {:.2} s\n", (i*100)/steps, n, (calc_t as f32)*0.001, (rend_t as f32)*0.001);
        }

    }
    Ok(())
}

