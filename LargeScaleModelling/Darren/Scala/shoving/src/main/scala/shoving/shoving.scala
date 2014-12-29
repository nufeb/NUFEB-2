/*
shoving.scala

*/

package shoving

import scala.swing._
import scala.swing.event._
import java.awt.{ Graphics2D, Color, BasicStroke }
import java.awt.geom.Ellipse2D
import java.awt.image.{ BufferedImage, WritableRaster }
import scala.annotation.tailrec

object shoving extends SimpleSwingApplication {

  def top = new MainFrame {
    title = "Shoving demo"

    val cell = Cell(0, 0, 0, 1, 0)
    val cell2 = Cell(1, 0, 0, 1, 0)
    val cellPop = List(cell, cell, cell2, cell2)

    val onButton = new Button {
      text = "Start"
    }
    val offButton = new Button {
      text = "Stop"
    }
    val panel = ImagePanel(1000, 800)
    drawImage(cellPop, panel.bi)
    contents = new BoxPanel(Orientation.Vertical) {
      contents += onButton
      contents += offButton
      contents += panel
      border = Swing.EmptyBorder(30, 30, 10, 30)
    }
    var cp = cellPop
    listenTo(onButton)
    listenTo(offButton)
    val timer = new javax.swing.Timer(1, Swing.ActionListener(e =>
      {
        cp = step(cp)
        cp = step(cp)
        drawImage(cp, panel.bi)
        println("%d cells".format(cp.length))
        panel.repaint()
      }))
    reactions += {
      case ButtonClicked(b) => {
        if (b.text == "Start")
          timer.start()
        else
          timer.stop()
      }
    }
  }

  @tailrec def evolve(cp: List[Cell], n: Int): List[Cell] = {
    println("%d cells".format(cp.length))
    //println(cp)
    if (n > 0)
      evolve(step(cp), n - 1)
    else
      cp
  }

  def step(cp: List[Cell]): List[Cell] = {
    val newCp = cp.par.map { _.drift }.map { _.grow }.map { _.age }.flatMap { _.divide }.flatMap { _.die }
    val shovedCp = newCp.map { c => c.shift(cp.map { _.force(c) }.reduce(_.add(_))) }
    shovedCp.toList
  }

  // Function to draw the cell pop on the BufferedImage
  def drawImage(cp: List[Cell], bi: BufferedImage): Unit = {
    val w = bi.getWidth
    val h = bi.getHeight
    val scaleFactor = 50 // height of window on original scale
    def scale(x: Double): Double = x * h / scaleFactor
    val big = bi.createGraphics()
    big.setColor(Color.white)
    big.fillRect(0, 0, w, h)
    big.setColor(Color.black)
    cp map { c => big.draw(new Ellipse2D.Double(h / 2 + scale(c.x), h / 2 + scale(c.y), scale(math.pow(c.s, 1.0 / 3)), scale(math.pow(c.s, 1.0 / 3)))) }
  }

  // ImagePanel class
  case class ImagePanel(bi: BufferedImage) extends Panel {
    override def paintComponent(g: Graphics2D) = {
      g.drawImage(bi, 0, 0, null)
    }
  }
  object ImagePanel {
    def apply(w: Int, h: Int) = {
      val bi = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB)
      val ip = new ImagePanel(bi)
      ip.preferredSize = new Dimension(w, h)
      ip
    }
  }

}


/* eof */


