package comp557.a2;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

import javax.swing.JPanel;
import javax.swing.border.TitledBorder;
import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;


import mintools.parameters.DoubleParameter;
import mintools.swing.VerticalFlowPanel;

/** 
 * Left Mouse Drag Arcball
 * @author kry
 */
/**** Marwen Kraiem 260955202 ****/

public class ArcBall {
		
	private DoubleParameter fit = new DoubleParameter( "Fit", 1, 0.5, 2 );
	private DoubleParameter gain = new DoubleParameter( "Gain", 1, 0.5, 2, true );
	
	private Vector3d v2 ;
	private Vector3d v1 ;


	/** The accumulated rotation of the arcball */
	Matrix4d R = new Matrix4d();

	public ArcBall() {
		R.setIdentity();
	}
	
	/** 
	 * Convert the x y position of the mouse event to a vector for your arcball computations 
	 * @param e
	 * @param v
	 */
	public Vector3d setVecFromMouseEvent( MouseEvent e ) {
		Component c = e.getComponent();
		Dimension dim = c.getSize();
		double width = dim.getWidth();
		double height = dim.getHeight();
		int mousex = e.getX();
		int mousey = e.getY();
		
		// TODO: Objective 1: finish arcball vector helper function
		double radius = Math.min(width, height)/this.fit.getFloatValue();
		double centerx = width/2.0d;
		double centery = height/2.0d;

		Vector3d v = new  Vector3d((mousex - centerx)/radius, -(mousey - centery)/radius,0.0d);
		double r = v.x*v.x + v.y*v.y;
	
		if (r > 1.0d){
			double s = 1.0/Math.sqrt(r);
			v.x = s*v.x;
			v.y = s*v.y;
			v.z = 0.0d;
		}
		else {v.z = Math.sqrt(1.0d-r);}

		return v;

	}
		
	public void attach( Component c ) {
		c.addMouseMotionListener( new MouseMotionListener() {		
			@Override
			public void mouseMoved( MouseEvent e ) {}
			@Override
			public void mouseDragged( MouseEvent e ) {				
				if ( (e.getModifiersEx() & MouseEvent.BUTTON1_DOWN_MASK) != 0 ) {
					// TODO: Objective 1: Finish arcball rotation update on mouse drag when button 1 down!
					v2 = setVecFromMouseEvent(e);	
					Vector3d axis = new Vector3d();
					axis.cross(v1, v2);
					if (axis.length()>1e-9) {	
						axis.normalize();
						AxisAngle4d aa = new AxisAngle4d( axis.x, axis.y, axis.z, gain.getDefaultValue()*v1.angle(v2) );
						//accumulated rotation
						Matrix4d tmpMatrix4d = new Matrix4d();
						tmpMatrix4d.set(aa);
						R.mul(tmpMatrix4d);
					}
				}
				v1.set(v2);
			}
		});
		c.addMouseListener( new MouseListener() {
			@Override
			public void mouseReleased( MouseEvent e) {}
			@Override
			public void mousePressed( MouseEvent e) {
				// TODO: Objective 1: arcball interaction starts when mouse is clicked
				v1 = setVecFromMouseEvent(e);
				
				
			}
			@Override
			public void mouseExited(MouseEvent e) {}
			@Override
			public void mouseEntered(MouseEvent e) {}
			@Override
			public void mouseClicked(MouseEvent e) {}
		});
	}
	
	
	public JPanel getControls() {
		VerticalFlowPanel vfp = new VerticalFlowPanel();
		vfp.setBorder( new TitledBorder("ArcBall Controls"));
		vfp.add( fit.getControls() );
		vfp.add( gain.getControls() );
		return vfp.getPanel();
	}
		
}
