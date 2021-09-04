//package varieties;
//
//import java.math.BigInteger;
//import java.util.Iterator;
//
//import fields.exceptions.InfinityException;
//import fields.helper.CoordinateRing;
//import fields.interfaces.Element;
//import fields.interfaces.Field;
//import fields.interfaces.Ideal;
//import fields.interfaces.Polynomial;
//import fields.polynomials.PolynomialIdeal;
//
//public class ProjectiveVarietyTmp<T extends Element<T>> implements VarietyTmp<T> {
//	private ProjectiveSpace<T> space;
//	private PolynomialIdeal<T> definingIdeal;
//	private CoordinateRing<T> coordinateRing;
//
//	public ProjectiveVarietyTmp(ProjectiveSpace<T> space, PolynomialIdeal<T> definingIdeal) {
//		this.space = space;
//		this.definingIdeal = definingIdeal;
//		this.coordinateRing = new CoordinateRing<T>(space.getRing(), this.definingIdeal);
//	}
//
//	@Override
//	public Exactness exactness() {
//		return space.exactness();
//	}
//
//	public final ProjectiveSpace<T> getSpace() {
//		return this.space;
//	}
//	
//	public final PolynomialIdeal<T>getIdeal() {
//		return definingIdeal;
//	}
//	
//	public final CoordinateRing<T> getHomogenousCoordinateRing() {
//		return coordinateRing;
//	}
//
//	@Override
//	public ProjectivePoint<T> getRandomElement() {
//		throw new UnsupportedOperationException();
//	}
//
//	@Override
//	public boolean isFinite() {
//		return this.space.isFinite() || this.coordinateRing.krullDimension() == 0;
//	}
//
//	@Override
//	public BigInteger getNumberOfElements() throws InfinityException {
//		throw new UnsupportedOperationException();
//	}
//
//	@Override
//	public Iterator<ProjectivePoint<T>> iterator() {
//		throw new UnsupportedOperationException();
//	}
//
//	@Override
//	public Field<T> getField() {
//		return this.space.getField();
//	}
//	
//	@Override
//	public int dimension() {
//		return coordinateRing.krullDimension() - 1;
//	}
//
//	@Override
//	public boolean hasRationalPoint(ProjectivePoint<T> p) {
//		return this.space.asIdeal(p).contains(this.definingIdeal);
//	}
//
//	public boolean contains(ProjectiveVarietyTmp<T> var) {
//		return var.definingIdeal.contains(this.definingIdeal);
//	}
//
//	public ProjectiveVarietyTmp<T> intersect(ProjectiveVarietyTmp<T> var) {
//		return new ProjectiveVarietyTmp<T>(this.space, this.space.getRing().add(this.definingIdeal, var.definingIdeal));
//	}
//
//	public ProjectiveVarietyTmp<T> unite(ProjectiveVarietyTmp<T> var) {
//		return new ProjectiveVarietyTmp<T>(this.space,
//				this.space.getRing().intersect(this.definingIdeal, var.definingIdeal));
//	}
//}
