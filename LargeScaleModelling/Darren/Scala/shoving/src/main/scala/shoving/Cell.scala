package shoving

import breeze.stats.distributions.{ Uniform, Gaussian }

case class Cell(x: Double, y: Double, z: Double, s: Double, a: Double) {

  def drift: Cell = {
    val dc = 0.01
    Cell(x + Gaussian(0, dc).draw, y + Gaussian(0, dc).draw, z + Gaussian(0, dc).draw, s, a)
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
    else {
      val jit = Gaussian(0, 0.1).sample(3)
      List(Cell(x - jit(0), y - jit(1), z - jit(2), 1.0, a / 4), Cell(x + jit(0), y + jit(1), z + jit(2), 1.0, 0))
    }
  }

  def die: List[Cell] = {
    if (Uniform(0, 1).draw < 0.0005 * math.log(a + 1))
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

  def rotate(th: Double): Cell = Cell(x * math.cos(th) - z * math.sin(th), y, x * math.sin(th) + z * math.cos(th), s, a)

}

case class Force(dx: Double, dy: Double, dz: Double) {
  def add(f: Force): Force = Force(dx + f.dx, dy + f.dy, dz + f.dz)
}
