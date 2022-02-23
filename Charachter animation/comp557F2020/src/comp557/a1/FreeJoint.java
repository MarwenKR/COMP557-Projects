package comp557.a1;

import com.jogamp.opengl.GLAutoDrawable;

import mintools.parameters.DoubleParameter;

//**** Marwen Kraiem 260955202 ****//

public class FreeJoint extends GraphNode {

	DoubleParameter tx;
	DoubleParameter ty;
	DoubleParameter tz;
	DoubleParameter rx;
	DoubleParameter ry;
	DoubleParameter rz;
		
	public FreeJoint( String name ) {
		super(name);
		dofs.add( tx = new DoubleParameter( name+" tx", 0, -3, 3 ) );		
		dofs.add( ty = new DoubleParameter( name+" ty", 0, -3, 3 ) );
		dofs.add( tz = new DoubleParameter( name+" tz", 0, -3, 3 ) );
		dofs.add( rx = new DoubleParameter( name+" rx", 0, -6.28, 6.28 ) );		
		dofs.add( ry = new DoubleParameter( name+" ry", 0, -6.28, 6.28 ) );
		dofs.add( rz = new DoubleParameter( name+" rz", 0, -6.28, 6.28 ) );
	}
	
	@Override
	public void display( GLAutoDrawable drawable, BasicPipeline pipeline ) {
		pipeline.push();
		
		// TODO: Objective 3: Freejoint, transformations must be applied before drawing children
		pipeline.translate(tx.getFloatValue(),ty.getFloatValue(), tz.getFloatValue());
		pipeline.rotate(rx.getFloatValue(), 1,0,0);
		pipeline.rotate(ry.getFloatValue(), 0,1,0);
		pipeline.rotate(rz.getFloatValue(), 0,0,1);
		
		super.display( drawable, pipeline );		
		pipeline.pop();
	}
	
}
