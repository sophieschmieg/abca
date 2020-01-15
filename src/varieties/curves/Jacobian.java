package varieties.curves;

import java.util.Map;
import java.util.TreeMap;

import varieties.curves.DivisorGroup.Divisor;
import varieties.curves.Jacobian.DivisorClass;
import fields.Element;
import fields.Group;
import fields.InfinityException;

public class Jacobian<T extends Element> implements Group<DivisorClass<T>> {
	private SmoothCurve<T> curve;
	private DivisorGroup<T> group;
	private Map<Divisor<T>, DivisorClass<T>> classes;
	
	public Jacobian(SmoothCurve<T> curve) {
		this.curve = curve;
		this.group = new DivisorGroup<T>();
		this.classes = new TreeMap<Divisor<T>, DivisorClass<T>>();
	}

	@Override
	public DivisorClass<T> getRandomElement() {
		return this.getDivisorClass(this.group.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return this.curve.isFinite();
	}

	@Override
	public int getNumberOfElements() throws InfinityException {
		return 0;
	}

	@Override
	public Iterable<DivisorClass<T>> getElements() throws InfinityException {
		return null;
	}

	@Override
	public DivisorClass<T> neutral() {
		return this.getDivisorClass(this.group.neutral());
	}

	@Override
	public DivisorClass<T> inverse(DivisorClass<T> t) {
		return this.getDivisorClass(this.group.inverse(t.getDivisor()));
	}

	@Override
	public DivisorClass<T> operate(DivisorClass<T> t1, DivisorClass<T> t2) {
		return this.getDivisorClass(this.group.operate(t1.getDivisor(), t2.getDivisor()));
	}
	public DivisorClass<T> getDivisorClass(Divisor<T> div) {
		if (this.classes.containsKey(div))
			return this.classes.get(div);
		DivisorClass<T> cl = new DivisorClass<T>(curve, div);
		for (DivisorClass<T> other : this.classes.values()) {
			if (other.equals(cl)) {
				this.classes.put(div, other);
				return this.getDivisorClass(div);
			}
		}
		this.classes.put(div, cl);
		return this.getDivisorClass(div);
	}
	public static class DivisorClass<T extends Element> implements Element {
		private SmoothCurve<T> curve;
		private Divisor<T> divisor;
		private DivisorGroup<T> group;
		
		private DivisorClass(SmoothCurve<T> curve, Divisor<T> divisor) {
			this.group = new DivisorGroup<T>();
			this.curve = curve;
			this.divisor = divisor;
	/*		if (this.divisor.getDegree() == 0 && this.curve instanceof EllipticCurve) {
				EllipticCurve<T> ell = (EllipticCurve<T>) this.curve;
				List<ProjectivePoint<T>> zeroes = new ArrayList<ProjectivePoint<T>>();
				List<ProjectivePoint<T>> poles = new ArrayList<ProjectivePoint<T>>();
				zeroes.addAll(this.divisor.getZeroes());
				poles.addAll(this.divisor.getPoles());
			}*/
			if (this.curve.isPrincipal(this.divisor)) {
				this.divisor = this.group.neutral();
			}
		}
		public Divisor<T> getDivisor() {
			return divisor;
		}
		@Override
		public boolean equals(Object o) {
			if (!(o instanceof DivisorClass))
				return false;
			@SuppressWarnings("unchecked")
			DivisorClass<T> dc = (DivisorClass<T>)o;
			if (this.curve.isPrincipal(this.group.operate(this.getDivisor(), this.group.inverse(dc.getDivisor()))))
				return true;
			return false;
		}
		@Override
		public String toString() {
			return this.divisor.toString();
		}
		@Override
		public int compareTo(Element o) {
			// TODO Auto-generated method stub
			return 0;
		}
	}
}
