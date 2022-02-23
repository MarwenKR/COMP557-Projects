package comp557.a1;

import javax.vecmath.Tuple3d;
import com.jogamp.opengl.GLAutoDrawable;
import mintools.parameters.DoubleParameter;

//**** Marwen Kraiem 260955202 ****//


public class RotaryJoint extends GraphNode {
	
	// angle
	DoubleParameter r;
	double minR;
	double maxR;

	// axis
	double ax;
	double ay;
	double az;
	
	// position
	double posx;
	double posy;
	double posz;
	
		
	public RotaryJoint( String name, double ax, double ay, double az, double posx, double posy, double posz, double minR, double maxR) {
		super(name);
		this.ax = ax;
		this.ay = ay;
		this.az = az;
		
		this.posx = posx;
		this.posy = posy;
		this.posz = posz;
		
		this.minR=minR;
		this.maxR=maxR;

		
		dofs.add( r = new DoubleParameter( name+" r", 0, minR, maxR ) );		
		
	}
		
	public void setPosition(Tuple3d t) {
		 posx = t.x ;
		 posy = t.y ;
		 posz = t.z ;
	}

	public void setAxis(Tuple3d t) {
		 ax = t.x ;
		 ay = t.y ;
		 az = t.z ;
	}
	
	@Override
	public void display( GLAutoDrawable drawable, BasicPipeline pipeline ) {
		pipeline.push();
		
		// TODO: Objective 4: RotaryJoint, transformations must be applied before drawing children
		pipeline.translate(posx,posy,posz);
		pipeline.rotate(r.getFloatValue(), ax, ay, az);
		
		pipeline.setModelingMatrixUniform(drawable.getGL().getGL4());
		super.display( drawable, pipeline );		
		pipeline.pop();
	}


}
