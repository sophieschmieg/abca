package varieties.curves;

import java.math.BigInteger;
import java.util.Iterator;
import java.util.Map;

import fields.exceptions.InfinityException;
import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.MathSet;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import varieties.RationalFunction;
import varieties.affine.AffineSubCover;
import varieties.curves.LineBundle.GlobalSection;
import varieties.projective.ProjectiveScheme;

public class LineBundle<T extends Element<T>> implements MathSet<GlobalSection<T>> {
	public static class GlobalSection<T extends Element<T>> extends AbstractElement<GlobalSection<T>> {
		private RationalFunction<T> sectionOnFirstCover;
		private LineBundle<T> lineBundle;

		private GlobalSection(LineBundle<T> lineBundle, RationalFunction<T> sectionOnFirstCover) {
			this.sectionOnFirstCover = sectionOnFirstCover;
			this.lineBundle = lineBundle;
		}

		@Override
		public int compareTo(GlobalSection<T> o) {
			return sectionOnFirstCover.compareTo(o.sectionOnFirstCover);
		}
	}

	private ProjectiveScheme<T> scheme;
	private AffineSubCover<T> cover;
	private Map<Integer, Map<Integer, CoordinateRingElement<T>>> transitionFunctions;

	public LineBundle(ProjectiveScheme<T> scheme, AffineSubCover<T> cover,
			Map<Integer, Map<Integer, CoordinateRingElement<T>>> transitionFunctions) {
		this.scheme = scheme;
		this.cover = cover;
		this.transitionFunctions = transitionFunctions;
	}

	@Override
	public Exactness exactness() {
		return getScheme().exactness();
	}

	public ProjectiveScheme<T> getScheme() {
		return scheme;
	}

	public AffineSubCover<T> getCover() {
		return cover;
	}

	@Override
	public boolean isFinite() {
		return scheme.getField().isFinite();
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new ArithmeticException("Not implemented");
	}

	@Override
	public Iterator<GlobalSection<T>> iterator() {
		throw new ArithmeticException("Not implemented");
	}

	@Override
	public GlobalSection<T> getRandomElement() {
		return new GlobalSection<>(this, scheme.getFunctionField().getRandomElement());
	}

	public CoordinateRingElement<T> transitionFunction(int firstCoverIndex, int secondCoverIndex) {
		return transitionFunctions.get(firstCoverIndex).get(secondCoverIndex);
	}

	public RationalFunction<T> onCover(GlobalSection<T> t, int coverIndex) {
		return scheme.getFunctionField().multiply(
				scheme.getFunctionField().getEmbedding(transitionFunction(0, coverIndex).getElement()),
				t.sectionOnFirstCover);
	}
}
