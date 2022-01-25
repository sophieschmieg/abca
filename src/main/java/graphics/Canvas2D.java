package graphics;

import java.awt.Dimension;
import java.awt.Graphics;
import java.util.ArrayList;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;

import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.integers.Integers;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.InnerProductSpace;
import fields.interfaces.Lattice;
import fields.interfaces.RealInnerProductSpace;
import fields.vectors.DualVectorSpace.Dual;
import fields.vectors.Polytope;
import fields.vectors.Vector;

public class Canvas2D<T extends Element<T>, S extends Element<S>> extends JPanel {
	private static final long serialVersionUID = 1L;
	private RealInnerProductSpace<T, S> space;
	private JFrame frame;
	private Reals r;
	private Real minX;
	private Real maxX;
	private Real minY;
	private Real maxY;
	private Polytope<T, S> asPolytope;

	private Lattice<?, T, S> lattice;
	private int minLatticeX;
	private int maxLatticeX;
	private int minLatticeY;
	private int maxLatticeY;

	private Polytope<T, S> polytope;

	public Canvas2D(RealInnerProductSpace<T, S> space, int sizeX, int sizeY, Real minX, Real maxX, Real minY,
			Real maxY) {
		this.space = space;
		if (space.dimension() != 2) {
			throw new ArithmeticException("Cannot draw non 2d things");
		}
		this.r = Reals.r(128);
		this.minX = minX;
		this.maxX = maxX;
		this.minY = minY;
		this.maxY = maxY;
		frame = new JFrame("drawing");
		frame.add(this);
		frame.setSize(new Dimension(sizeX, sizeY));
		frame.setVisible(true);
		InnerProductSpace<T, Dual<T, S>> dualSpace = space.getDual();
		List<Dual<T, S>> constraints = new ArrayList<>();
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(1)));
		constraints.add(dualSpace.negative(dualSpace.getUnitVector(2)));
		constraints.add(dualSpace.getUnitVector(1));
		constraints.add(dualSpace.getUnitVector(2));
		List<T> values = new ArrayList<>();
		values.add(space.fromReal(r.negative(minX)));
		values.add(space.fromReal(r.negative(minY)));
		values.add(space.fromReal(maxX));
		values.add(space.fromReal(maxY));
		asPolytope = new Polytope<>(space, constraints, values);
	}

	private int coordinateX(S point) {
		Vector<T> asVector = space.asVector(point);
		T t = asVector.get(1);
		int width = this.getWidth();
		Real adjusted = r.multiply(width, r.divide(r.subtract(space.asReal(t), minX), r.subtract(maxX, minX)));
		return adjusted.round().intValueExact();
	}

	private int coordinateY(S point) {
		Vector<T> asVector = space.asVector(point);
		T t = asVector.get(2);
		int height = this.getHeight();
		Real adjusted = r.multiply(height, r.divide(r.subtract(space.asReal(t), minY), r.subtract(maxY, minY)));
		return height - adjusted.round().intValueExact();
	}

	private boolean inCanvas(S point) {
		return asPolytope.contains(point);
	}

	public <R extends Element<R>> void setLattice(Lattice<R, T, S> lattice, int minX, int maxX, int minY, int maxY) {
		this.lattice = lattice;
		this.minLatticeX = minX;
		this.maxLatticeX = maxX;
		this.minLatticeY = minY;
		this.maxLatticeY = maxY;
		repaint();
	}

	public void setPolytope(Polytope<T, S> polytope) {
		this.polytope = polytope;
		repaint();
	}

	@Override
	public void paint(Graphics g) {
		if (this.lattice != null) {
			drawLattice(g, lattice, minLatticeX, maxLatticeX, minLatticeY, maxLatticeY);
		}
		if (this.polytope != null) {
			drawPolytope(g, polytope);
		}
	}

	private void drawPoint(Graphics graphics, S point) {
		if (!inCanvas(point)) {
			return;
		}
		graphics.fillRect(coordinateX(point) - 1, coordinateY(point) - 1, 3, 3);
	}

	private <R extends Element<R>> void drawLattice(Graphics graphics, Lattice<R, T, S> lattice, int minX, int maxX,
			int minY, int maxY) {
		Integers z = Integers.z();
		for (int i = minX; i <= maxX; i++) {
			for (int j = minY; j <= maxY; j++) {
				drawPoint(graphics,
						lattice.embedding(lattice.fromVector(new Vector<>(z.getInteger(i), z.getInteger(j)))));
			}
		}
	}

	private void drawConstraint(Graphics graphics, Dual<T, S> constraint, T value) {
		Vector<T> asVector = space.asVector(constraint.dual());
		Field<T> f = space.getField();
		List<S> possiblePoints = new ArrayList<>();
		if (!asVector.get(2).equals(f.zero())) {
			possiblePoints.add(space.fromVector(new Vector<>(space.fromReal(minX),
					f.divide(f.subtract(value, f.multiply(space.fromReal(minX), asVector.get(1))), asVector.get(2)))));
			possiblePoints.add(space.fromVector(new Vector<>(space.fromReal(maxX),
					f.divide(f.subtract(value, f.multiply(space.fromReal(maxX), asVector.get(1))), asVector.get(2)))));
		}
		if (!asVector.get(1).equals(f.zero())) {
			possiblePoints.add(space.fromVector(new Vector<>(
					f.divide(f.subtract(value, f.multiply(space.fromReal(minY), asVector.get(2))), asVector.get(1)),
					space.fromReal(minY))));
			possiblePoints.add(space.fromVector(new Vector<>(
					f.divide(f.subtract(value, f.multiply(space.fromReal(maxY), asVector.get(2))), asVector.get(1)),
					space.fromReal(maxY))));
		}
		List<S> points = new ArrayList<>();
		for (S point : possiblePoints) {
			if (inCanvas(point)) {
				points.add(point);
			}
		}
		if (points.size() != 2) {
			return;
		}
		graphics.drawLine(coordinateX(points.get(0)), coordinateY(points.get(0)), coordinateX(points.get(1)),
				coordinateY(points.get(1)));
	}

	private void drawPolytope(Graphics graphics, Polytope<T, S> polytope) {
		for (int i = 0; i < polytope.getConstraints().size(); i++) {
			drawConstraint(graphics, polytope.getConstraints().get(i), polytope.getRightHandSide().get(i));
		}
	}
}
