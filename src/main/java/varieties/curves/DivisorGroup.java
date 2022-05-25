package varieties.curves;

import java.math.BigInteger;
import java.util.Collections;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.integers.Integers;
import fields.interfaces.Element;
import fields.interfaces.Group;
import varieties.projective.ProjectivePoint;
import varieties.projective.ProjectiveScheme;

public class DivisorGroup<T extends Element<T>> implements Group<WeilDivisor<T>> {
	private ProjectiveScheme<T> scheme;

	DivisorGroup(ProjectiveScheme<T> scheme) {
		this.scheme = scheme;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public WeilDivisor<T> getRandomElement() {
		Integers z = Integers.z();
		int degree = z.getRandomElement(z.getInteger(11)).intValueExact() - 5;
		int number = z.getRandomElement(z.getInteger(10)).intValueExact();
		Map<ProjectivePoint<T>, Integer> result = new TreeMap<>();
		for (int i = 0; i < number + degree; i++) {
			ProjectivePoint<T> zero = scheme.getRandomElement();
			if (result.containsKey(zero)) {
				result.put(zero, result.get(zero) + 1);
			} else {
				result.put(zero, 1);
			}
		}
		for (int i = 0; i < number - degree; i++) {
			ProjectivePoint<T> pole = scheme.getRandomElement();
			if (result.containsKey(pole)) {
				result.put(pole, result.get(pole) - 1);
			} else {
				result.put(pole, -1);
			}
		}
		return new WeilDivisor<>(result);
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<WeilDivisor<T>> iterator() {
		throw new InfinityException();
	}

	@Override
	public WeilDivisor<T> neutral() {
		Map<ProjectivePoint<T>, Integer> map = Collections.emptyMap();
		return new WeilDivisor<>(map);
	}

	@Override
	public WeilDivisor<T> inverse(WeilDivisor<T> t) {
		Map<ProjectivePoint<T>, Integer> map = new TreeMap<ProjectivePoint<T>, Integer>();
		for (ProjectivePoint<T> p : t.getSupport())
			map.put(p, -t.get(p));
		return new WeilDivisor<>(map);
	}

	@Override
	public WeilDivisor<T> operate(WeilDivisor<T> t1, WeilDivisor<T> t2) {
		Map<ProjectivePoint<T>, Integer> map = new TreeMap<>();
		Set<ProjectivePoint<T>> keyset = new TreeSet<>();
		keyset.addAll(t1.getSupport());
		keyset.addAll(t2.getSupport());
		for (ProjectivePoint<T> p : keyset) {
			map.put(p, t1.get(p) + t2.get(p));
		}
		return new WeilDivisor<>(map);
	}
}
