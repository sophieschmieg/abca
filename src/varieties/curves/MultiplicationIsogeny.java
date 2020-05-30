package varieties.curves;

import java.math.BigInteger;

import fields.Element;
import varieties.ProjectivePoint;

public class MultiplicationIsogeny<T extends Element> implements Isogeny<T> {
	private BigInteger n;
	private EllipticCurve<T> curve;
	
	public MultiplicationIsogeny(BigInteger n, EllipticCurve<T> curve) {
		this.n = n;
		this.curve = curve;
	}

	@Override
	public EllipticCurve<T> getDomain() {
		return curve;
	}

	@Override
	public EllipticCurve<T> getRange() {
		return curve;
	}

	@Override
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
		return curve.multiply(n, point);
	}

	@Override
	public BigInteger getDegree() {
		return n.multiply(n);
	}

	@Override
	public Isogeny<T> getDual() {
		return this;
	}

}
