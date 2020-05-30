package varieties.curves;

import java.math.BigInteger;

import fields.Element;
import varieties.ProjectivePoint;

public class CompositionIsogeny<T extends Element> implements Isogeny<T> {
	private Isogeny<T> firstIsogeny;
	private Isogeny<T> secondIsogeny;

	public CompositionIsogeny(Isogeny<T> firstIsogeny, Isogeny<T> secondIsogeny) {
		if (!firstIsogeny.getRange().equals(secondIsogeny.getDomain())) {
			throw new ArithmeticException("Isogenies cannot be composed!");
		}
		this.firstIsogeny = firstIsogeny;
		this.secondIsogeny = secondIsogeny;
	}

	@Override
	public EllipticCurve<T> getDomain() {
		return firstIsogeny.getDomain();
	}

	@Override
	public EllipticCurve<T> getRange() {
		return secondIsogeny.getRange();
	}

	@Override
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
		return secondIsogeny.evaluate(firstIsogeny.evaluate(point));
	}

	@Override
	public BigInteger getDegree() {
		return firstIsogeny.getDegree().multiply(secondIsogeny.getDegree());
	}

	@Override
	public Isogeny<T> getDual() {
		return new CompositionIsogeny<T>(secondIsogeny.getDual(), firstIsogeny.getDual());
	}

}
