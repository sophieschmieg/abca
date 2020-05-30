package varieties.curves;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import fields.Element;
import fields.Field;
import varieties.ProjectivePoint;

public class Frobenius<T extends Element> implements Isogeny<T> {
	private EllipticCurve<T> domain;
	private EllipticCurve<T> range;
	private Field<T> field;
	private BigInteger p;
	
	public Frobenius(EllipticCurve<T> domain) {
		this.domain = domain;
		this.field = domain.getField();
		this.p = this.field.characteristic();
		this.range = new EllipticCurve<>(field, field.power(domain.getA(), p), field.power(domain.getB(), p));
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
			coords.add(field.power(coord, p));
		}
		return new ProjectivePoint<T>(field, coords);
	}

	@Override
	public BigInteger getDegree() {
		return p;
	}

	@Override
	public Isogeny<T> getDual() {
		return new Isogeny<T>() {

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
					Set<T> roots = field.roots(coord, p.intValueExact());
					if (roots.size() != 1) {
						throw new ArithmeticException("pth root not unique");
					}
					coords.addAll(roots);
				}
				return domain.multiply(p, new ProjectivePoint<T>(field, coords));
			}

			@Override
			public BigInteger getDegree() {
				return p;
			}

			@Override
			public Isogeny<T> getDual() {
				return Frobenius.this;
			}
		};
	}

}
