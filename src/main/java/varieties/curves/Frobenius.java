package varieties.curves;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class Frobenius<T extends Element<T>> implements Isogeny<T> {
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
			PolynomialRing<T> polynomials = domain.asProjectiveVariety().homogenousPolynomialRing();
			for (int var = 0; var < polynomials.numberOfVariables(); var++) {
				asPolynomials.add(polynomials.getVarPower(var + 1, q.intValueExact()));
			}
			asMorphism = new ProjectiveMorphism<>(domain.asProjectiveVariety(), range.asProjectiveVariety(), asPolynomials);
		}
		return asMorphism;
	}

	@Override
	public Isogeny<T> getDual() {
		return new Isogeny<T>() {
			private ProjectiveMorphism<T> asMorphism = null;

			@Override
			public EllipticCurve<T> getDomain() {
				return range;
			}

			@Override
			public EllipticCurve<T> getRange() {
				return domain;
			}

			@Override
			public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
				List<T> coords = new ArrayList<>();
				for (T coord : point.getCoords()) {
					Map<T, Integer> roots = field.roots(coord, q.intValueExact());
					if (roots.size() != 1) {
						throw new ArithmeticException("qth root not unique");
					}
					coords.addAll(roots.keySet());
				}
				return domain.multiply(q, new ProjectivePoint<T>(field, coords));
			}

			@Override
			public BigInteger getDegree() {
				return q;
			}

			@Override
			public Isogeny<T> getDual() {
				return Frobenius.this;
			}

			@Override
			public ProjectiveMorphism<T> asMorphism() {
				if (asMorphism == null) {
					List<Polynomial<T>> asPolynomials = new ArrayList<>();
					Isogeny<T> mulQ = new MultiplicationIsogeny<T>(q, range);
					for (Polynomial<T> mulQPolynomial : mulQ.asMorphism().asPolynomials()) {
						asPolynomials.add(range.asProjectiveVariety().homogenousPolynomialRing().characteristicRoot(mulQPolynomial, power));
					}
					asMorphism = new ProjectiveMorphism<>(range.asProjectiveVariety(), domain.asProjectiveVariety(), asPolynomials);
				}
				return asMorphism;
			}
		};
	}
}
