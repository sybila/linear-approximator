import com.github.daemontus.Formula
import com.github.daemontus.ModelFile
import com.github.daemontus.tokenizer.Rules
import com.github.daemontus.tokenizer.parseLines
import java.io.File
import java.util.*

fun List<*>.printArray() = this.joinToString(prefix = "(", postfix = ")", separator = ", ")

fun Double.safeString(): String {
    return String.format(Locale.ROOT, "%f", this)
}

fun main(args: Array<String>) {
    println("Arguments: function_file point_count segment_count fast")
    val file = File(args[0])
    val modelFile = file.useLines {
        it.toList().parseLines()
    }
    val vars = modelFile.declarations.filterIsInstance<ModelFile.Declaration.Variable>()
    if (vars.size != 1) error("More than one variable!")
    val variable = vars[0]
    val functions = modelFile.declarations.filterIsInstance<ModelFile.Declaration.Function>()
    if (functions.any { it.arguments.isNotEmpty() }) error("Functions can't have arguments!")
    val bounds = variable.bounds.first.toEvaluable(null)(0.0) to variable.bounds.second.toEvaluable(null)(0.0)
    val evals = functions.map { it.value.toEvaluable(variable.name) }

    val thresholds = computeThresholds(bounds, args[1].toInt(), args[2].toInt(), evals, args[3].toBoolean())
    println("${variable.name}: ${thresholds.toList().printArray()}")
    for ((f, e) in functions.zip(evals)) {
        println("${f.name}: ${thresholds.zip(thresholds.map(e)).map { (a, b) -> "[${a.safeString()}, ${b.safeString()}]" }.printArray()}")
    }
}

typealias Evaluable = (Double) -> Double

fun Formula.toEvaluable(varName: String?): Evaluable {
    return when (this) {
        is Formula.Number -> { _ -> this.value }
        is Formula.Text -> error("Text not allowed!")
        is Formula.Function -> {
            val args = this.args.map { it.toEvaluable(varName) }
            when (this.name) {
                Rules.Math.Addition.id -> { x: Double -> args.fold(0.0) { a, b -> a + b(x) } }
                Rules.Math.Multiplication.id -> { x: Double -> args.fold(1.0) { a, b -> a * b(x) } }
                Rules.Math.Subtraction.id -> {
                    if (args.size == 1) {
                        { x: Double -> -1 * args[0](x) }
                    } else {
                        { x: Double -> args.drop(1).fold(args[0](x)) { a, b -> a - b(x) } }
                    }
                }
                Rules.Math.Division.id -> {
                    if (args.size == 1) {
                        { x: Double -> 1 / args[0](x) }
                    } else {
                        { x: Double -> args.drop(1).fold(args[0](x)) { a, b -> a / b(x) } }
                    }
                }
                "mod" -> {
                    val f: Evaluable = if (args.size == 2) {
                        { x: Double ->
                            var a = args[0](x)
                            val b = args[1](x)
                            if (b < 0 || a < 0) error("Negative values not supported in mod")
                            while (a > b) a -= b
                            a
                        }
                    } else error("Modulo operation has two arguments")
                    f
                }
                "tanh" -> {
                    val f: Evaluable = if (args.size == 1) {
                        { x: Double ->
                            Math.tanh(args[0](x))
                        }
                    } else error("tanh has one argument")
                    f
                }
                "pow" -> {
                    val f: Evaluable = if (args.size == 2) {
                        { x: Double ->
                            Math.pow(args[0](x), args[1](x))
                        }
                    } else error("pow has two arguments")
                    f
                }
                varName -> { x: Double -> x }
                else -> error("Unsupported function ${this.name}")
            }
        }
    }
}

//compute approximated thresholds for one variable
private fun computeThresholds(bounds: Pair<Double, Double>, pointCount: Int, segmentCount: Int, functions: List<Evaluable>, fast: Boolean): DoubleArray {

    val xPoints = findEvaluationPoints(pointCount, bounds)
    val curves = Array(functions.size) { f -> DoubleArray(pointCount) { functions[f](xPoints[it]) } }

    return if (fast) {
        fastLinearApproximation(xPoints, curves, segmentCount)
    } else {
        linearApproximation(xPoints, curves, segmentCount)
    }

}

//compute points at which the function should be evaluated
private fun findEvaluationPoints(pointCount: Int, bounds: Pair<Double, Double>): DoubleArray {
    val (min, max) = bounds
    val dx = (max - min) / (pointCount-1)
    return DoubleArray(pointCount) { i -> min + dx * i }
}

//performs approximation using accurate cost function
private fun linearApproximation(xPoints: DoubleArray, curves: Array<DoubleArray>, segmentCount: Int): DoubleArray {

    val hCost = Array(xPoints.size) { n -> DoubleArray(xPoints.size) { i ->
        if (i > n-2 || i == 0) {
            0.0 //no one will read this!
        } else {
            curves.maxByDoubleIndexed { ic -> segmentError(xPoints, curves[ic], i, n) }
        }
    }}

    return approximation(segmentCount, xPoints.size, xPoints, curves, hCost)
}

//performs approximation using simplified cost function
private fun fastLinearApproximation(xPoints: DoubleArray, curves: Array<DoubleArray>, segmentCount: Int): DoubleArray {

    val sy = Array(curves.size) { DoubleArray(xPoints.size) { 0.0 } }
    val sy2 = Array(curves.size) { DoubleArray(xPoints.size) { 0.0 } }
    val sxy = Array(curves.size) { DoubleArray(xPoints.size) { 0.0 } }

    val sx = DoubleArray(xPoints.size) { 0.0 }
    val sx2 = DoubleArray(xPoints.size) { 0.0 }

    for (ic in curves.indices) {
        sy2[ic][0] = curves[ic][0] * curves[ic][0]
        sy [ic][0] = curves[ic][0]
        sxy[ic][0] = curves[ic][0] * xPoints[0]

        for (ip in 1 until xPoints.size) {
            sy2[ic][ip] = sy2[ic][ip-1] + (curves[ic][ip] * curves[ic][ip])
            sy [ic][ip] = sy [ic][ip-1] + (curves[ic][ip])
            sxy[ic][ip] = sxy[ic][ip-1] + (curves[ic][ip] * xPoints[ip])
        }
    }

    sx2[0] = xPoints[0] * xPoints[0]
    sx[0] = xPoints[0]

    for (ip in 1 until xPoints.size) {
        sx2[ip] = sx2[ip-1] + (xPoints[ip] * xPoints[ip])
        sx [ip] = sx [ip-1]  + xPoints[ip]
    }

    val hCost = Array(xPoints.size) { n -> DoubleArray(xPoints.size) { i ->
        if (i > n-2 || i == 0) {
            0.0 //no one will read this!
        } else {
            curves.maxByDoubleIndexed { ic ->
                val a = (curves[ic][n] - curves[ic][0]) / (xPoints[n] - xPoints[0])
                val b = (curves[ic][0] * xPoints[n] - curves[ic][n] * xPoints[0]) / (xPoints[n] - xPoints[0])
                (sy2[ic][n] - sy2[ic][i-1]) - 2 * a * (sxy[ic][n] - sxy[ic][i-1]) - 2 * b * (sy[ic][n] - sy[ic][i-1]) + a * a * (sx2[n] - sx2[i-1]) + 2 * a * b * (sx[n] - sx[i-1]) + b * (n - i)
            }
        }
    }}

    return approximation(segmentCount, xPoints.size, xPoints, curves, hCost)
}

//performs approximation process using given cost function
private fun approximation(
        segmentCount: Int,
        pointCount: Int,
        points: DoubleArray,
        values: Array<DoubleArray>,
        costs: Array<DoubleArray>
): DoubleArray {

    val father = Array(segmentCount) { IntArray(pointCount) { 0 } }

    val cost = DoubleArray(pointCount - 1) { i ->
        values.maxByDouble { curve ->
            segmentError(points, curve, 0, i + 1)
        }
    }

    for (m in 1 until segmentCount) {

        for (n in (pointCount - 1).downTo(2)) {

            var minError = cost[n-2]
            var minIndex = n - 1

            for (i in m..(n-2)) {

                val currentError = cost[i-1] + costs[n][i]

                if (currentError < minError) {
                    minError = currentError
                    minIndex = i
                }
            }

            cost[n-1] = minError
            father[m][n] = minIndex

        }

    }


    val results = DoubleArray(segmentCount+1) { 0.0 }

    var pointIndex = pointCount - 1
    results[segmentCount] = points[pointIndex]

    for (i in (segmentCount-1).downTo(0)) {
        pointIndex = father[i][pointIndex]
        results[i] = points[pointIndex]
    }

    return results
}

private fun segmentError(x: DoubleArray, y: DoubleArray, first: Int, last: Int): Double {
    // Compute line segment coefficients
    val a = (y[last] - y[first]) / (x[last] - x[first])
    val b = (y[first] * x[last] - y[last] * x[first]) / (x[last] - x[first])

    // Compute error for the line segment
    var e = 0.0
    @Suppress("LoopToCallChain")    // much faster this way
    for (k in first..last) {
        e += (y[k] - a * x[k] - b) * (y[k] - a * x[k] - b)
    }
    e /= (a*a + 1)

    return e
}

//Utility functions used when computing the approximated model
private inline fun <T> Array<T>.maxByDouble(action: (T) -> Double): Double {
    var max = Double.NEGATIVE_INFINITY
    @Suppress("LoopToCallChain")    // much faster this way
    for (e in this) {
        max = Math.max(max, action(e))
    }
    return max
}

private inline fun <T> Array<T>.maxByDoubleIndexed(action: (Int) -> Double): Double {
    var max = Double.NEGATIVE_INFINITY
    @Suppress("LoopToCallChain")    // much faster this way
    for (e in this.indices) {
        max = Math.max(max, action(e))
    }
    return max
}
