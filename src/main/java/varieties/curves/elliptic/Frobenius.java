package varieties.curves.elliptic;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class Frobenius<T extends Element<T>> extends AbstractElement<Isogeny<T>> implements Isogeny<T> {
	private EllipticCurve<T> domain;
	private EllipticCurve<T> range;
	private Field<T> field;
	private BigInteger q;
	private int power;
	private ProjectiveMorphism<T> asMorphism;

	public Frobenius(EllipticCurve<T> domain, int power) {
		this.domain = domain;
		this.field = domain.getField();
		this.q = this.field.characteristic().pow(power);
		this.power = power;
		this.range = new EllipticCurve<>(field, field.power(domain.getA(), q), field.power(domain.getB(), q));
	}

	@Override
	public String toString() {
		return "F_" + q;
	}

	@Override
	public int compareTo(Isogeny<T> o) {
		if (o instanceof Frobenius<?>) {
			Frobenius<T> other = (Frobenius<T>) o;
			int cmp = domain.compareTo(other.domain);
			if (cmp != 0) {
				return cmp;
			}
			return power - other.power;
		}
		return getClass().getName().compareTo(o.getClass().getName());
	}

	@Override
	public EllipticCurve<T> getDomain() {
		return domain;
	}

	@Override
	public EllipticCurve<T> getRange() {
		return range;
	}

	@Override
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
		List<T> coords = new ArrayList<>();
		for (T coord : point.getCoords()) {
			coords.add(field.power(coord, q));
		}
		return new ProjectivePoint<T>(field, coords);
	}

	@Override
	public BigInteger getDegree() {
		return q;
	}

	@Override
	public ProjectiveMorphism<T> asMorphism() {
		if (asMorphism == null) {
			List<Polynomial<T>> asPolynomials = new ArrayList<>();
			PolynomialRing<T> polynomials = domain.asGenericProjectiveScheme().homogenousPolynomialRing();
			for (int var = 0; var < polynomials.numberOfVariables(); var++) {
				asPolynomials.add(polynomials.getVarPower(var + 1, q.intValueExact()));
			}
			asMorphism = new ProjectiveMorphism<>(domain.asGenericProjectiveScheme(), range.asGenericProjectiveScheme(),
					asPolynomials);
		}
		return asMorphism;
	}

	private static class DualFrobenius<T extends Element<T>> extends AbstractElement<Isogeny<T>> implements Isogeny<T> {
		private ProjectiveMorphism<T> asMorphism = null;
		private Frobenius<T> frobenius;

		private DualFrobenius(Frobenius<T> frobenius) {
			this.frobenius = frobenius;
		}

		@Override
		public String toString() {
			return frobenius.toString() + "~";
		}

		@Override
		public int compareTo(Isogeny<T> o) {
			if (o instanceof DualFrobenius<?>) {
				DualFrobenius<T> other = (DualFrobenius<T>) o;
				return frobenius.compareTo(other.frobenius);
			}
			return getClass().getName().compareTo(o.getClass().getName());
		}

		@Override
		public EllipticCurve<T> getDomain() {
			return frobenius.range;
		}

		@Override
		public EllipticCurve<T> getRange() {
			return frobenius.domain;
		}

		@Override
		public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
			List<T> coords = new ArrayList<>();
			for (T coord : point.getCoords()) {
				Map<T, Integer> roots = frobenius.field.roots(coord, frobenius.q.intValueExact());
				if (roots.size() != 1) {
					throw new ArithmeticException("qth root not unique");
				}
				coords.addAll(roots.keySet());
			}
			return frobenius.domain.multiply(frobenius.q, new ProjectivePoint<T>(frobenius.field, coords));
		}

		@Override
		public BigInteger getDegree() {
			return frobenius.q;
		}

		@Override
		public Isogeny<T> getDual() {
			return frobenius;
		}

		@Override
		public ProjectiveMorphism<T> asMorphism() {
			if (asMorphism == null) {
				List<Polynomial<T>> asPolynomials = new ArrayList<>();
				Isogeny<T> mulQ = new MultiplicationIsogeny<T>(frobenius.q, frobenius.range);
				for (Polynomial<T> mulQPolynomial : mulQ.asMorphism().asPolynomials()) {
					asPolynomials.add(frobenius.range.asGenericProjectiveScheme().homogenousPolynomialRing()
							.characteristicRoot(mulQPolynomial, frobenius.power));
				}
				asMorphism = new ProjectiveMorphism<>(frobenius.range.asGenericProjectiveScheme(),
						frobenius.domain.asGenericProjectiveScheme(), asPolynomials);
			}
			return asMorphism;
		}

//		@Override
		public List<ProjectivePoint<T>> kernelGenerators() {
			if (frobenius.domain.isSupersingular()) {
				return Collections.emptyList();
			}
			return Collections.singletonList(frobenius.range.getRandomTorsionPoint(frobenius.q));
		}
	}

	@Override
	public Isogeny<T> getDual() {
		return new DualFrobenius<>(this);
	}

//	@Override
	public List<ProjectivePoint<T>> kernelGenerators() {
		return Collections.emptyList();
	}
}
