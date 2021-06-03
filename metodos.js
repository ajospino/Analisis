const math = require("mathjs/lib/browser/math");
const ndarray = require('ndarray');
var zeros = require('zeros'),
    fill = require('ndarray-fill');
 

class methods{

  static biseccion(fx, a, b, maxNum) {
    console.log("Biseccion: ");
    let cond = 0.0000001;
    let i = 1;
    let xm = 0;
    let xm0 = 0;
    let fa = 0;
    let fm = 0;
    let check = 0;
    let error = 1;

    console.log(
      "| i  |     a         |    xm         |      b        |     f(x)   |    E    |"
    );

    while (error > cond && i < maxNum) {
      //Primera iteracion

      if (i == 1) {
        xm = (a + b) / 2;
        fa = math.evaluate(fx);
        fm = fa(xm);
        fa = fa(a);
      } else {
        //Condicion de a y b
        if (check < 0) {
          b = xm;
        } else {
          a = xm;
        }

        xm0 = xm;
        xm = (a + b) / 2;
        fa = math.evaluate(fx);
        fm = fa(xm);
        fa = fa(a);
        error = math.abs(xm - xm0);
      }
      //Valor de check
      if (fm * fa < 0) {
        check = -1;
      } else {
        check = 1;
      }
      console.log(
        "|  " +
          i +
          " | " +
          a.toFixed(10) +
          "  | " +
          xm.toFixed(10) +
          "  | " +
          b.toFixed(10) +
          "  | " +
          fm.toPrecision(6) +
          "  | " +
          error +
          "  |"
      );
      i = i + 1;
    }
  }

  static reglaFalsa(fx, a, b, maxNum) {
    console.log("Regla Falsa:");
    let cond = 0.0000001;
    let i = 1;
    let error = 1;
    let xm = 0;
    let xm0 = 0;
    let f = 0;
    let fa = 0;
    let fb = 0;
    let fm = 0;

    console.log(
      "| i  |     a         |    xm         |     b         |      f(x)    |    E    |"
    );

    while (error > cond && i < maxNum) {
      if (i == 1) {
        f = math.evaluate(fx);
        fa = f(a);
        fb = f(b);
        xm = (fb * a - fa * b) / (fb - fa);
        fm = f(xm);
      } else {
        if (fa * fm < 0) {
          b = xm;
        } else {
          a = xm;
        }
        xm0 = xm;
        fa = f(a);
        fb = f(b);
        xm = (fb * a - fa * b) / (fb - fa);
        fm = f(xm);
        error = math.abs(xm - xm0);
      }
      console.log(
        "|  " +
          i +
          " | " +
          a.toFixed(10) +
          "  | " +
          xm.toFixed(10) +
          "  | " +
          b.toFixed(10) +
          "  | " +
          fm.toPrecision(6) +
          "  | " +
          error +
          "  |"
      );
      i = i + 1;
    }
  }

  static puntoFijo(fx, gx, x0, maxNum) {
    console.log("Punto Fijo:");
    let cond = 0.0000001;
    let error = 1;
    let i = 0;
    let x = x0;
    let f = 0;
    let g = 0;

    console.log(
      "| i  |          xi          |       g(xi)          |    f(xi)  |     E   |"
    );

    while (error > cond && i < maxNum) {
      if (i == 0) {
        f = math.evaluate(fx);
        f = f(x);
        g = math.evaluate(gx);
        g = g(x);
        console.log(
          "|  " +
            i +
            " | " +
            x.toFixed(16) +
            "  | " +
            g +
            "  | " +
            f.toPrecision(6) +
            "  | " +
            error.toPrecision(6) +
            " |    "
        );
      } else {
        x = g;
        f = math.evaluate(fx);
        f = f(x);
        g = math.evaluate(gx);
        g = g(x);
        error = math.abs(x - x0);
        x0 = x;
        console.log(
          "|  " +
            i +
            " | " +
            x +
            "  | " +
            g +
            "  | " +
            f.toPrecision(6) +
            "  | " +
            error.toPrecision(6) +
            "  |    "
        );
      }

      i = i + 1;
    }
  }

  static secante(x0, x1, fx, maxNum) {
    console.log("Secante:");
    let error = 1;
    let i = 0;
    let cond = 0.0000001;
    let f0 = 0;
    let f1 = 0;
    let x = x0;

    console.log(
      "| iter |         xi          |          f(x)       |          E          |"
    );

    while (error > cond && i < maxNum) {
      if (i == 0) {
        f0 = math.evaluate(fx);
        f0 = f0(x0);
        console.log(
          "|  " +
            i +
            "   | " +
            x.toFixed(16) +
            "  | " +
            f1.toFixed(16) +
            "  | " +
            error.toFixed(16) +
            "  |"
        );
      } else if (i == 1) {
        f1 = math.evaluate(fx);
        f1 = f1(x1);
        console.log(
          "|  " +
            i +
            "   | " +
            x1.toFixed(16) +
            "  | " +
            f1 +
            "  | " +
            error.toFixed(16) +
            "  |"
        );
      } else {
        x = x1;
        x1 = x1 - (f1 * (x1 - x0)) / (f1 - f0);
        x0 = x;
        f0 = f1;
        f1 = math.evaluate(fx);
        f1 = f1(x1);
        error = math.abs(x1 - x0);
        console.log(
          "|  " + i + "   | " + x1 + "  | " + f1 + "  | " + error + "  |"
        );
      }

      i = i + 1;
    }
  }
/** 
  static gaussSimple(Ma, b) {
    const matrixA = math.matrix(Ma);
    const vectorB = b;

    const n = matrixA.size()[0];

    let M;
    math.evaluate(`M = [matrixA, vectorB]`, {
      M,
      matrixA,
      vectorB,
    });

    for (let i = 0; i < n - 1; i++) {
      for (let j = i + 1; j < n; j++) {
        if (math.subset(M, math.index(i, j)) !== 0) {
          math.evaluate(`M[j,i:n+1]=M[j,i:n+1]-((M[j,i]/M[i,i])*M[i,i:n+1])`, {
            M,
            n,
            i,
            j,
          });
        }
      }
    }
  }
  */
  
  static vanderMon(x,y) {
        
  if( x.dimension !== 1 ) {
    throw new TypeError('vandermonde: error: x debe ser un vector de una dimensión');
  }
  let inicio= new Date();

  setTimeout(function(){
    terminarScript(inicio);
  },0);

  var i;
  var M = x.shape[0];
  
  var N = arguments[2] || M;
  var reversed = !! arguments[3];
  var dtype = arguments[4] || 'float64';
  var v = zeros([M,N],dtype);
  if( ! reversed ) {
    fill(v,function(i,j) {
      return Math.pow( x.get(i), j );
    });
  } else {
    fill(v,function(i,j) {
      return Math.pow( x.get(M-1-i), j);
    });
  }
  
  var vres = ndarray(new Float64Array(v.size),[math.sqrt(v.size),math.sqrt(v.size)]);
  var cont = x.size-1;
  var cont2 = 0;
    for(let u = 0;u < v.size; u++){
      if(cont< 0){
        cont = x.size-1;
      }
      
      if(u%x.size == 0 && u!=0){
        
        cont2 = cont2+x.size;
      }
      // console.log(v.data[cont2+cont])
      
      vres.data[u] = v.data[cont2+cont];
      
      cont--;
    }
    
    const textMatrix = math.ones(x.size,x.size);
    var matrixres = ""
    
    for(let p=0;p<vres.size/x.size;p++){
      for(let q=0;q<vres.size/x.size;q++){
        textMatrix.subset(math.index(q,p),vres.get(q,p))
       
       matrixres = matrixres + " " + vres.get(p,q) 
      }
      matrixres = matrixres + "\n"
    }
    
    var coeffRes = "";
    var coeff = (math.multiply(math.inv(textMatrix),y))
    var cont1 = x.size-1;
    for(let ty = 0; ty<coeff._size[0];ty++){
      if(ty == coeff._size[0]-1){
        coeffRes = coeffRes + coeff._data[ty] + "x^" + cont1 
      }
      else{
        coeffRes = coeffRes + coeff._data[ty] + "x^" + cont1+ " + " 
      }
    
    cont1--;
    }
    console.log("La Matriz de vandermonde es:")
    console.log()
    console.log(matrixres)
    console.log()
    console.log("El polinomio de la matriz es : ")
    console.log(coeffRes)
    console.log()
    function terminarScript(fechaInicial){
      let fin=new Date();
      let diferencia = fin-fechaInicial;
      console.log("El tiempo que se demora es: " + diferencia + " milisegundos")
    }
    return vres;
        
        
}
  
    //No funcional
    static difDivid2(Vx,Vy) {
      n= math.size(Vx)
      var X = Vx
    
      D=math.zeros(n[0],n[0])
    
      Y = math.transpose(Vy)
      for(let i = 0; i<n; i++ ){
    
        D._data[0,i][0] = Y[i]
        console.log()
      }
    
      for(let i = 1 ; i<n[0];i++){
    
        var aux0 = D._data[math.range(i-1,n[0]),i-1]
        aux = math.diff(aux0)
    
        aux2 = math.subtract(math.subset(X,math.index(math.range(i,n[0]))), math.subset(X,math.index(math.range(0,n[0]-1))))
    
        // D._data[math.range(i,n[0]),i] = math.divide(aux,math.transpose(aux2))
      }
      console.log(aux)
      console.log(math.transpose(aux2))
    
      Coef = math.diag(D)
      // console.log(D._data)
      // console.log(Coef._data)
    }
 
    /**
  static difDivid(Vx,Vy){
    const n = Vx.length;
    const X = math.resize(math.matrix(Vx),[1,4]);
    const Y = math.matrix(Vy);
    let D = math.zeros(n,n);
    const YT = math.resize(Y,[4,1]);
    
    
    D = math.evaluate(`D[:,0]=YT`, {
        D,
        YT,
      });

    
    let aux0;
    for (let i = 0; i < n - 1; i++) {
      math.evaluate(`aux0 = D[i-1:n,i-1]`, {
        aux0,
        D,
        n,
        i,
      });
      aux = math.eval(math.derivative(aux0))
      math.evaluate(`aux2 = X[i:n] - X[0:n-i]`, {
        aux2,
        X,
        n,
        i,
      });
      math.evaluate(`D[i:n,i] = aux/math.transponse(aux2)`, {
        D,
        aux2,
        n,
        i,
      });
    }
    
  }
  */
  static steffenSen( x0,f,tol,nmax ) {
    let x1 = 0
    let x2 = 0
    var x= x0
    x = eval(f)
    x1 = x
    x2 = eval(f)
    var xn = x0 -math.pow((x1-x0),2)/(x2-2*x1+x0)
    let error =math.abs(x2-xn)
    let i = 1
    while(error>=tol && i < nmax){
            x0 = xn
            x =x0
            x1= eval(f)
            x=x1
            x2=eval(f)
            xn = x0 -math.pow((x1-x0),2)/(x2-2*x1+x0)
            error = math.abs(x0-xn);
            i++;
            console.log("Error en la iteracion " + i + " es de " +error)
            console.log("El valor de xn en la iteracion " + i + " es de " +xn)
            console.log()
    }
    console.log(x + " " + x1 )
  }
}

/** 
var k =([[4,-1,0,3],[1,15.5,3,8],[0,-1.3,-4,1.1],[14,5,-2,30]])
var ksp = math.sparse([[4,-1,0,3],[1,15.5,3,8],[0,-1.3,-4,1.1],[14,5,-2,30]])
math.lusolve(k,[1,1,1,1])
console.log(math.lusolve(k,[1,1,1,1]))
*/

// console.log(math.slu(ksp,0,0.0000001))
// console.log(math.lup(k))


//methods.reglaFalsa("f(x)=log(sin(x)^2 + 1) -(0.5)", 0, 1, 100)
//methods.secante(0.5, 1, "f(x) = log(sin(x)^2 + 1) -(0.5)", 100)
//methods.biseccion("f(x) = log(sin(x)^2 + 1) -(0.5)",0,1,100)
//methods.puntoFijo('f(x) = log(sin(x)^2 + 1) - (0.5) -x ', 'f(x) = log(sin(x)^2 + 1) -(0.5)', -0.5, 100)
/** 
var x = ndarray(new Float64Array([-1, 0, 3, 4]));
var h = [15.5, 3, 8, 1];
methods.vanderMon(x,h);
*/

//methods.steffenSen(1.5,'math.sqrt(10/(x+4))',0.00001,10)


//methods.difDivid2([-1,0,3,4],[15.5,3,8,1])