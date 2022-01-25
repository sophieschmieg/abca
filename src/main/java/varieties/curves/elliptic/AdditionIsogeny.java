package varieties.curves.elliptic;

import java.math.BigInteger;

import fields.helper.AbstractElement;
import fields.interfaces.Element;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class AdditionIsogeny<T extends Element<T>> extends AbstractElement<Isogeny<T>> implements Isogeny<T> {
	private EllipticCurve<T> domain;
	private EllipticCurve<T> range;
	private Isogeny<T> firstSummand;
	private Isogeny<T> secondSummand;
	private BigInteger degree;

	public AdditionIsogeny(Isogeny<T> firstSummand, Isogeny<T> secondSummand) {
		if (!firstSummand.getDomain().equals(secondSummand.getDomain())) {
			throw new ArithmeticException("domain differs");
		}
		if (!firstSummand.getRange().equals(secondSummand.getRange())) {
			throw new ArithmeticException("range differs");
		}
		this.domain = firstSummand.getDomain();
		this.range = firstSummand.getRange();
		this.firstSummand = firstSummand;
		this.secondSummand = secondSummand;
	}

	@Override
	public int compareTo(Isogeny<T> o) {
		if (o instanceof AdditionIsogeny<?>) {
			AdditionIsogeny<T> other = (AdditionIsogeny<T>) o;
			int cmp = firstSummand.compareTo(other.firstSummand);
			if (cmp != 0) {
				return cmp;
			}
			return secondSummand.compareTo(other.secondSummand);
		}
		return getClass().getName().compareTo(o.getClass().getName());
	}

	@Override
	public BigInteger getDegree() {
		if (degree == null) {
			CompositionIsogeny<T> composition = new CompositionIsogeny<>(getDual(), this);
			degree = domain.asMultiplicationIsogenyModOrder(composition);
		}
		return degree;
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
		return range.add(firstSummand.evaluate(point), secondSummand.evaluate(point));
	}

	@Override
	public ProjectiveMorphism<T> asMorphism() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Isogeny<T> getDual() {
		return new AdditionIsogeny<>(firstSummand.getDual(), secondSummand.getDual());
	}

}
