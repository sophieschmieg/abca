package varieties.curves;

import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.polynomials.CoordinateRing;
import varieties.FunctionField;
import varieties.RationalFunction;
import varieties.Scheme;
import varieties.projective.GenericProjectiveScheme;
import varieties.projective.ProjectivePoint;
import varieties.projective.ProjectiveScheme;

public interface SmoothCurve<T extends Element<T>> extends Scheme<T, ProjectivePoint<T>>, ProjectiveScheme<T> {
	public int getEmbeddingDimension();
	public GenericProjectiveScheme<T> asGenericProjectiveScheme();
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p);
	public List<Polynomial<T>> getTangentSpace(ProjectivePoint<T> p);
	public List<RationalFunction<T>> getRiemannRochSpace(WeilDivisor<T> div);
	public boolean isPrincipal(WeilDivisor<T> div);
	public CoordinateRing<T> getCoordinateRing();
	public FunctionField<T> getFunctionField();
}
