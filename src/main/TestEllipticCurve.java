package main;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.ArrayList;

import varieties.ProjectivePoint;
import varieties.curves.DivisorGroup;
import varieties.curves.DivisorGroup.Divisor;
import varieties.curves.EllipticCurve;
import fields.Element;
import fields.Field;
import fields.InfinityException;
import fields.RationalFunction;
import fields.Polynomial;

public class TestEllipticCurve<T extends Element> {
	private Field<T> field;
	private EllipticCurve<T> curve;
	private ProjectivePoint<T> x1;
	private ProjectivePoint<T> x2;
	private ProjectivePoint<T> r;
	private List<ProjectivePoint<T>> empty;
	private DivisorGroup<T> group;
	private Map<ProjectivePoint<T>, T> tau;
	private Map<ProjectivePoint<T>, Map<ProjectivePoint<T>, T>> sigma;
	
	/**
	 * @param args
	 * @throws InfinityException 
	 */
	@SuppressWarnings("unchecked")
	public TestEllipticCurve(Field<T> field) throws InfinityException {
		this.field = field;
		this.curve = new EllipticCurve<T>(field, field.one(), field.one());
                int n = (int)Math.ceil(4 * (double)Math.sqrt(field.getNumberOfElements()));
                for (int i = 0; i < 10; i++) {
                  System.out.println(i + " Division Polynomial: " + this.curve.getDivisionPolynomial(i));
                }
		int num = curve.getNumberOfElements();
                System.out.println("Number of elements: " + num);
		x1 = new ProjectivePoint<T>(field, this.field.getInteger(2), this.field.getInteger(1), field.one());
		x2 = new ProjectivePoint<T>(field, this.field.getInteger(2), this.field.getInteger(4), field.one());
		empty = Collections.emptyList();
		group = new DivisorGroup<T>();
		System.out.println(curve.getThirdIntersection(x2, x2));
		for (ProjectivePoint<T> p : curve.getElements()) {
			if (!p.equals(curve.neutral()) && 
					!p.equals(x1) &&
					!p.equals(x2) &&
					!p.equals(curve.add(x1, x1)) &&
					!p.equals(curve.add(x2, x2)) &&
					!p.equals(curve.add(x1, x2)) &&
					!x1.equals(curve.add(x2, p)) &&
					!x2.equals(curve.add(x1, p)) &&
					!curve.neutral().equals(curve.add(x1, p)) &&
					!curve.neutral().equals(curve.add(x2, p))) {
				r = p;
				break;
			}
		}
		System.out.println("Field: " + field);
		System.out.println("Field Characteristic: " + field.characteristic());
		System.out.println("Number of Field Elements: " + field.getNumberOfElements());
		System.out.println("Curve: " + curve);
		System.out.println("Glueing points: " + x1 + " " + x2);
		System.out.println("Alternative Divisorpole: " + r);
                List<Integer> div = new ArrayList<>();
                for (int i = 1; i <= num; i++) {
                  if (num % i == 0) {
                    div.add(i);
                  }
                }
                int counter = 1;
		for (ProjectivePoint<T> p : curve.getElements()) {
			int i = -1;
			for (int d : div) {
			  if (curve.multiply(d, p).equals(curve.neutral())) {
                            i = d;
                            break;
                          }
                        }
                        if (i > 2) {
                          Polynomial<T> divPoly = curve.getDivisionPolynomial(i);
                          System.out.println(i + " division polynomial evaluated at " + p + ": " + divPoly.evaluate(p.getCoords()));
                        }
//			T t = tau(p);
			System.out.println(counter + " Point: " + p + " Order: " + i); //+ " Tau: " + t/ + " Root: " + this.field.hasRoot(t, i));
                        counter++;
/*			if (!this.field.hasRoot(t, i))
				System.err.println("Points not working!");*/
		}
		/*
		for (ProjectivePoint<T> p : curve.getElements()) {
			for (ProjectivePoint<T> q : curve.getElements()) {
				System.out.println(p + " + " + q + " = " + curve.add(p, q) + " " + sigma(p, q) + " " + field.multiply(field.multiply(tau(p), tau(q)), field.inverse(tau(curve.add(p, q)))));
			}
		}*/
	}
	private T tau(ProjectivePoint<T> p) {
		if (this.tau == null)
			this.tau = new TreeMap<ProjectivePoint<T>, T>();
		if (this.tau.containsKey(p))
			return this.tau.get(p);
		ProjectivePoint<T> q = p;
		T res = field.one();
		//int i = 1;
		while (!q.equals(this.curve.neutral())) {
			res = field.multiply(res, sigma(p, q));
			q = curve.add(p, q);
			//i++;
		}
		this.tau.put(p, res);
		return res; //TODO: ite Wurzel! 
	}
	private T sigma(ProjectivePoint<T> p, ProjectivePoint<T> q) {
		if (this.sigma == null)
			this.sigma = new TreeMap<ProjectivePoint<T>, Map<ProjectivePoint<T>,T>>();
		if (this.sigma.containsKey(p) && this.sigma.get(p).containsKey(q))
			this.sigma.get(p).get(q);
		Divisor<T> d1;
		Divisor<T> d2;
		Divisor<T> d3;
		if (p.equals(x1))
			d1 = new Divisor<T>(Collections.singletonList(curve.add(x1, r)), Collections.singletonList(r));
		else if (p.equals(x2))
			d1 = new Divisor<T>(Collections.singletonList(curve.add(x2, r)), Collections.singletonList(r));
		else
			d1 = new Divisor<T>(Collections.singletonList(p), Collections.singletonList(curve.neutral()));
		if (q.equals(x1))
			d2 = new Divisor<T>(Collections.singletonList(curve.add(x1, r)), Collections.singletonList(r));
		else if (q.equals(x2))
			d2 = new Divisor<T>(Collections.singletonList(curve.add(x2, r)), Collections.singletonList(r));
		else
			d2 = new Divisor<T>(Collections.singletonList(q), Collections.singletonList(curve.neutral()));
		if (curve.add(p, q).equals(x1) || curve.add(p, q).equals(x2))
			d3 = new Divisor<T>(Collections.singletonList(r), empty);
		else
			d3 = new Divisor<T>(Collections.singletonList(curve.neutral()), empty);
		Divisor<T> div = group.operate(group.operate(d1, d2), d3);
		RationalFunction<T> f = curve.getRiemannRochSpace(div).get(0);
		T res = field.multiply(f.evaluate(x1).getDehomogenisedCoord(1, 2), f.evaluate(x2).getDehomogenisedCoord(2, 1));
		if (!this.sigma.containsKey(p))
			this.sigma.put(p, new TreeMap<ProjectivePoint<T>, T>());
		this.sigma.get(p).put(q, res);
		return res;
	}

}
