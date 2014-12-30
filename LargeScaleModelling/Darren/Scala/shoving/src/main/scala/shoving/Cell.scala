package shoving

import breeze.stats.distributions.Uniform

case class Cell(x: Double, y: Double, z: Double, s: Double, a: Double) {

  def drift: Cell = {
    Cell(x + Uniform(-0.001, 0.001).draw, y + Uniform(-0.001, 0.001).draw, z + Uniform(-0.001, 0.001).draw, s, a)
  }

  def grow: Cell = {
    Cell(x, y, z, s + Uniform(0.001, 0.01).draw, a)
  }

  def age: Cell = {
    Cell(x, y, z, s, a + 1)
  }

  def divide: List[Cell] = {
    if (s < 2.0)
      List(this)
    else
      // TODO: randomly displace
      List(Cell(x - 0.1, y, z, 1.0, a/2), Cell(x + 0.1, y, z, 1.0, 0))
  }

  def die: List[Cell] = {
    if (Uniform(0, 1).draw < 0.0001 * math.log(a + 1))
      List()
    else
      List(this)
  }

  def force(c: Cell): Force = {
    val dx = x - c.x
    val dy = y - c.y
    val dz = z - c.z
    val dr = math.sqrt(dx * dx + dy * dy + dz * dz)
    val td = math.pow(s, 1.0 / 3) + math.pow(c.s, 1.0 / 3)
    if (dr > 2 * td)
      Force(0, 0, 0)
    else {
      //val f = 0.01 * (dr - td) * (dr - 2 * td)
      //val f = if (dr < td) 0.02 else -0.001
      val f = 0.01 * (if (dr < td) (10.0 - 5.0 * dr / td) else (-2.0 + dr / td))
      Force(f * dx / dr, f * dy / dr, f * dz / dr)
    }
  }

  def shift(f: Force): Cell = Cell(x - f.dx, y - f.dy, z - f.dz, s, a)

  def rotate(th: Double): Cell = Cell(x*math.cos(th)-z*math.sin(th),y,x*math.sin(th)+z*math.cos(th),s,a)
  
}

case class Force(dx: Double, dy: Double, dz: Double) {
  def add(f: Force): Force = Force(dx + f.dx, dy + f.dy, dz + f.dz)
}
