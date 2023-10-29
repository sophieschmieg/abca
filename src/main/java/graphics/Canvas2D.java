package graphics;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
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
	private class Point {
		private Point(S point, Color color) {
			this.point = point;
			this.color = color;
		}

		private S point;
		private Color color;
	}

	private class Arrow {
		private Arrow(S origin, S vector, Color color) {
			this.origin = origin;
			this.vector = vector;
			this.color = color;
		}

		private S origin;
		private S vector;
		private Color color;
	}

	private static final long serialVersionUID = 1L;
	private RealInnerProductSpace<T, S> space;
	private JFrame frame;
	private Object lock;
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

	private List<Point> markedPoints;
	private List<Arrow> markedVectors;
	private boolean coordinates;

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
		lock = new Object();
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
		markedPoints = new ArrayList<>();
		markedVectors = new ArrayList<>();
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

	public void setCoordinates(boolean coordinates) {
		this.coordinates = coordinates;
		repaint();
	}

	public <R extends Element<R>> void setLattice(Lattice<R, T, S> lattice, int minX, int maxX, int minY, int maxY) {
		this.lattice = lattice;
		this.minLatticeX = minX;
		this.maxLatticeX = maxX;
		this.minLatticeY = minY;
		this.maxLatticeY = maxY;
		lattice.getModuleGenerators();
		repaint();
	}

	public void setPolytope(Polytope<T, S> polytope) {
		this.polytope = polytope;
		repaint();
	}

	public void addMarkedPoint(S point, Color color) {
		this.markedPoints.add(new Point(point, color));
		repaint();
	}

	public void addMarkedVector(S start, S vector, Color color) {
		this.markedVectors.add(new Arrow(start, vector, color));
		repaint();
	}

	public void waitUntilClose() {
		Thread t = new Thread() {
			@Override
			public void run() {
				synchronized (lock) {
					while (frame.isVisible())
						try {
							lock.wait();
						} catch (InterruptedException e) {
							e.printStackTrace();
						}
				}
			}
		};
		t.start();
		frame.addWindowListener(new WindowAdapter() {

			@Override
			public void windowClosing(WindowEvent arg0) {
				synchronized (lock) {
					frame.setVisible(false);
					lock.notify();
				}
			}

		});
		try {
			t.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void paint(Graphics g) {
		drawCoordinates(g);
		if (this.lattice != null) {
			drawLattice(g, lattice, minLatticeX, maxLatticeX, minLatticeY, maxLatticeY);
		}
		if (this.polytope != null) {
			drawPolytope(g, polytope);
		}
		for (Arrow vector : markedVectors) {
			Color c = g.getColor();
			g.setColor(vector.color);
			drawArrow(g, vector.origin, vector.vector);
			g.setColor(c);
		}
		for (Point point : markedPoints) {
			Color c = g.getColor();
			g.setColor(point.color);
			drawPoint(g, point.point);
			g.setColor(c);
		}
	}

	private void drawPoint(Graphics graphics, S point) {
		if (!inCanvas(point)) {
			return;
		}
		graphics.fillRect(coordinateX(point) - 1, coordinateY(point) - 1, 3, 3);
	}

	private void drawArrow(Graphics graphics, S start, S vector) {
		if (!inCanvas(start)) {
			return;
		}
		S dest = space.add(start, vector);
		if (!inCanvas(dest)) {
			return;
		}
		graphics.drawLine(coordinateX(start), coordinateY(start), coordinateX(dest), coordinateY(dest));
		int angle = 15;
		Real arrowAngle = r.multiply(r.pi(), r.getDouble(angle / 180.0));
		int length = 10;
		Real scaleX = r.multiply(length, r.divide(r.subtract(maxX, minX), r.getInteger(this.getWidth())));
		Real scaleY = r.multiply(length, r.divide(r.subtract(maxY, minY), r.getInteger(this.getHeight())));
		Vector<T> asVector = space.asVector(vector);
		Real vectorAngle = r.add(r.pi(), r.arctan2(space.asReal(asVector.get(2)), space.asReal(asVector.get(1))));
		Vector<Real> arrow1 = r.cosineAndSine(r.add(vectorAngle, arrowAngle));
		S point1 = space.add(dest, space.fromVector(new Vector<>(space.fromReal(r.multiply(scaleX, arrow1.get(1))),
				space.fromReal(r.multiply(scaleY, arrow1.get(2))))));
		Vector<Real> arrow2 = r.cosineAndSine(r.subtract(vectorAngle, arrowAngle));
		S point2 = space.add(dest, space.fromVector(new Vector<>(space.fromReal(r.multiply(scaleX, arrow2.get(1))),
				space.fromReal(r.multiply(scaleY, arrow2.get(2))))));
		graphics.drawLine(coordinateX(point1), coordinateY(point1), coordinateX(dest), coordinateY(dest));
		graphics.drawLine(coordinateX(point2), coordinateY(point2), coordinateX(dest), coordinateY(dest));
	}

	private void drawCoordinates(Graphics graphics) {
		if (coordinates) {
			Real epsX = r.divide(r.subtract(maxX, minX), r.getInteger(getWidth()));
			Real epsY = r.divide(r.subtract(maxY, minY), r.getInteger(getHeight()));
			S minXPoint = space
					.fromVector(new Vector<>(space.fromReal(r.add(minX, epsX)), space.getValueField().zero()));
			S maxXPoint = space
					.fromVector(new Vector<>(space.fromReal(r.subtract(maxX, epsX)), space.getValueField().zero()));
			drawArrow(graphics, minXPoint, space.subtract(maxXPoint, minXPoint));
			S minYPoint = space
					.fromVector(new Vector<>(space.getValueField().zero(), space.fromReal(r.add(minY, epsY))));
			S maxYPoint = space
					.fromVector(new Vector<>(space.getValueField().zero(), space.fromReal(r.subtract(maxY, epsY))));
			drawArrow(graphics, minYPoint, space.subtract(maxYPoint, minYPoint));
			FontMetrics metrics = graphics.getFontMetrics();
			S unit1 = space.getUnitVector(1);
			drawPoint(graphics, unit1);
			graphics.drawString("1", coordinateX(unit1) - metrics.stringWidth("1") / 2,
					coordinateY(unit1) + metrics.getHeight());
			S unit2 = space.getUnitVector(2);
			drawPoint(graphics, unit2);
			graphics.drawString("1", coordinateX(unit2) - metrics.stringWidth("1 "),
					coordinateY(unit2) + metrics.getHeight() / 3);
		}
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
