package fields.quaternions;

import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;

public class RealQuaternions extends AbstractQuaternions<Real> {
	private Reals base;
	
	private RealQuaternions(Reals base, Real a, Real b) {
		super(base, a, b);
		this.base = base;
	}

	public static RealQuaternions quaternions(Reals r) {
		return new RealQuaternions(r, r.getInteger(-1), r.getInteger(-1));
	}

	public static RealQuaternions quaternions(Reals r, Real a, Real b) {
		return new RealQuaternions(r, a, b);
	}

	@Override
	public int hilbertSymbol() {
		if (a().compareTo(base.zero()) < 0 && b().compareTo(base.zero()) < 0) {
			return -1;
		}
		return 1;
	}

	@Override
	public Real discriminant() {
		return hilbertSymbol() == -1 ? base.getInteger(2) : base.one();
	}

	@Override
	public boolean isIntegral() {
		return hilbertSymbol() == -1;
	}

	@Override
	protected Quaternion<Real> splittingElement() {
		if (a().compareTo(base.zero()) > 0) {
			return i();
		} else if (b().compareTo(base.zero()) > 0) {
			return j();
		}
		throw new ArithmeticException("Not a split quaternion algebra!");
	}

	@Override
	protected Quaternion<Real> normalize(Quaternion<Real> t) {
		Real norm = reducedNorm(t);
		if (norm.equals(base.zero())) {
			return t;
		}
		Real sqrt = base.positiveSqrt(norm);
		return divide(t, getEmbedding(sqrt));
	}

}
