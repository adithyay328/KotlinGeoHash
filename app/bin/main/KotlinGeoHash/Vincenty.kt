package KotlinGeoHash

// Some constants used by Vincenty's Formulae(https://en.wikipedia.org/wiki/Vincenty%27s_formulae)
val a = 6378137.0
val f = 1/298.257223563
val b = (1-f) * a

fun vincentyDist(phiOne : Double, L1 : Double, phiTwo : Double, L2: Double, maxIters : Int = 100, minUpdate : Double = Math.pow(10.0, -12.0)) : Triple<Double, Double, Double> {
  // Some values that are constants w.r.t to the paramaters. More performant to just compute
  // them once

  // Want L as a double in radians
  val L = Math.toRadians(L2 - L1)

  // Also make sure to cast phiOne and 2 into radians
  val U1 = Math.atan( (1 - f) * Math.tan(Math.toRadians(phiOne) ))
  val U2 = Math.atan( (1 - f) * Math.tan(Math.toRadians(phiTwo) ))

  val cosU1 = Math.cos(U1)
  val sinU1 = Math.sin(U1)
  
  val cosU2 = Math.cos(U2)
  val sinU2 = Math.sin(U2)
  
  // Lambda represents distance; initialize to L
  var lambda = L

  // Current number of iterations
  var currIters = 0
  
  // Loop only breaks when sufficient precision is met or
  // max number of iterations have run
  while(true) {
    val sinLambda = Math.sin(lambda)
    val cosLambda = Math.cos(lambda)

    val sinSigma = Math.sqrt(  Math.pow(cosU2 * sinLambda, 2.0) + Math.pow( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda , 2.0) )

    val cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda

    val sigma = Math.atan(sinSigma / cosSigma)

    val sinAlpha = ( cosU1 * cosU2 * sinLambda ) / sinSigma

    val sinSquaredAlpha = Math.pow(sinAlpha, 2.0)

    val cosSquaredAlpha = 1 - sinSquaredAlpha
    val cosTwoSigmaM = cosSigma - ( (2 * sinU1 * sinU2) / cosSquaredAlpha )

    val cosSquaredTwoSigmaM = Math.pow(cosTwoSigmaM, 2.0)
    val C = (f / 16) * cosSquaredAlpha * (4 + f * (4 - (3 * cosSquaredAlpha) ))
    val newLambda = L + (1 - C ) * f * sinAlpha * ( sigma + C * sinSigma * (cosTwoSigmaM + C * cosSigma * (-1 + 2 * cosSquaredTwoSigmaM)) )

    currIters++

    // How much lambda has changed since the last update.
    // Once this update is small enough, we've converged. Keep
    // iterating until we get there or we run out of iterations
    val lambdaUpdate = Math.abs(newLambda - lambda)
    lambda = newLambda

    // Final calculations
    if((currIters > maxIters) or (lambdaUpdate <= minUpdate)) {
      val uSquared = cosSquaredAlpha * ( (Math.pow(a, 2.0) - Math.pow(b, 2.0)) / Math.pow(b, 2.0))
      
      val A = 1 + (uSquared / 16384) * (4096 + uSquared * (-768 + uSquared * (320 - 175 * uSquared)))
      
      val B = ( uSquared / 1024 ) * ( 256 + uSquared * (-128 + uSquared * (74 - 47 * uSquared)) )

      val deltaSigma = B * sinSigma * (cosTwoSigmaM + (B/4) * ( cosSigma * (-1 + 2 * cosSquaredTwoSigmaM) ) - (B / 6) * cosTwoSigmaM * (-3 + 4 * Math.pow(sinSigma, 2.0)) * ( -3 + 4 * cosSquaredTwoSigmaM ) )
      
      val s = b * A * (sigma - deltaSigma)
      
      val alpha1 = Math.atan2(cosU2 * sinLambda, cosU1 * sinU2 - sinU1 * cosU2 * cosLambda)
      
      val alpha2 = Math.atan2(cosU1 * sinLambda, -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda)

      return Triple(s, alpha1, alpha2)
    }

  }
}