package comp557.a1;

import javax.vecmath.Tuple2d;
import javax.vecmath.Tuple3d;

import com.jogamp.opengl.GLAutoDrawable;

import mintools.parameters.DoubleParameter;

//**** Marwen Kraiem 260955202 ****//


public class SphericalJoint extends GraphNode {

	// Euler angles
	DoubleParameter rx;
	DoubleParameter ry;
	DoubleParameter rz;
	double minRx, maxRx;
	double minRy, maxRy;
	double minRz, maxRz;
	
	// position
	double posx;
	double posy;
	double posz;
		
	public SphericalJoint(  String name, double posx, double posy, double posz, double minRx, double maxRx, double minRy, double maxRy, double minRz, double maxRz) {
		super(name);
		
		
		this.posx = posx;
		this.posy = posy;
		this.posz = posz;
		
		this.minRx=minRx;
		this.maxRx=maxRx;
		this.minRy=minRy;
		this.maxRy=maxRy;
		this.minRz=minRz;
		this.maxRz=maxRz;
		
		dofs.add( rx = new DoubleParameter( name+" rx", (minRx+maxRx)/2, minRx, maxRx ) );		
		dofs.add( ry = new DoubleParameter( name+" ry", (minRy+maxRy)/2, minRy, maxRy ) );
		dofs.add( rz = new DoubleParameter( name+" rz", (minRz+maxRz)/2, minRz, maxRz ) );
	}
	
	public void setPosition(Tuple3d t) {
		 posx = t.x ;
		 posy = t.y ;
		 posz = t.z ;
	}
	
	public void setMinMaxRx(Tuple2d t) {
		 minRx = t.x ;
		 maxRz = t.y ;
	}
	
	public void setMinMaxRy(Tuple2d t) {
		 minRy = t.x ;
		 maxRy = t.y ;
	}
	
	public void setMinMaxRz(Tuple2d t) {
		 minRz = t.x ;
		 maxRz = t.y ;
	}


	@Override
	public void display( GLAutoDrawable drawable, BasicPipeline pipeline ) {
		pipeline.push();
		
		pipeline.translate(posx,posy,posz);
		
		// TODO: Objective 5: SphericalJoint, transformations must be applied before drawing children
		pipeline.rotate(rx.getFloatValue(), 1,0,0);
		pipeline.rotate(ry.getFloatValue(), 0,1,0);
		pipeline.rotate(rz.getFloatValue(), 0,0,1);
		
		super.display( drawable, pipeline );		
		pipeline.pop();
	}


	


	
}
