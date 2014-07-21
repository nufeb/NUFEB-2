import scala.swing._
import scala.swing.event._
import java.awt.{Graphics2D,Color,BasicStroke}
import scala.math._
import scala.util.Random
import scala.annotation.tailrec

object MySwingApp extends SimpleSwingApplication {

 def top = new MainFrame {
  title="Swing app"
  val button=new Button {
   text="Generate"
  }
  val panel=new Canvas {
   preferredSize=new Dimension(500,500)
  }
  contents=new BoxPanel(Orientation.Vertical) {
   contents+=button
   contents+=panel
   border=Swing.EmptyBorder(30,30,10,30)
  }
  listenTo(button)
  var nClicks=0
  reactions+={
   case ButtonClicked(b) =>
    panel.repaint()
  }
 }


}


class Canvas extends Panel {

 override def paintComponent(g: Graphics2D) = {
  g.clearRect(0,0,size.width,size.height)
  g.setColor(Color.blue)
  // g.plotPoint(400,300)
  val r=Random
  agglom(size.width/2,0,100,r,g)
 }

 @tailrec final def agglom(x0: Int,y0: Int,n: Int,r: Random,g: Graphics2D): Unit = {
  if (n>0) {
    g.drawLine((r.nextDouble*100).toInt,(r.nextDouble*100).toInt,(r.nextDouble*100).toInt,(r.nextDouble*100).toInt)
    agglom(x0,y0,n-1,r,g)
  }
 }



}




