
rooms.scene2 = function() {

lib3D2();

description = `<b>Scene 2</b>
               <p>
               Hierarchic 3D scene
               <br>
               with triangle meshes.
	       <p> <input type=range id=arm_length value=28> arm length
	       <br><input type=range id=leg_length value=40> leg length
	       <br><input type=range id=rate> rate
	       `;

code = {
'init':`





















   S.redPlastic    = [.2,.1,.1,0,  .5,.2,.2,0,  2,2,2,20,  0,0,0,0];
   S.greenPlastic  = [.1,.2,.1,0,  .2,.5,.2,0,  2,2,2,20,  0,0,0,0];
   S.bluePlastic   = [.1,.1,.2,0,  .2,.2,.5,0,  2,2,2,20,  0,0,0,0];
   S.whitePlastic  = [.2,.2,.2,0,  .5,.5,.5,0,  2,2,2,20,  0,0,0,0];
   S.purplePlastic = [.5,.2,.5,0,  .2,.1,.2,0,  2,2,2,20,  0,0,0,0];
   S.orangePlastic = [.5,.5,.1,0,  .2,.2,0,0,   2,2,2,20,  0,0,0,0];



   // A SQUARE IS A TRIANGLE MESH WITH JUST TWO TRIANGLES

   S.squareMesh = [ -1, 1, 0,  0,0,1,  0,1,
                     1, 1, 0,  0,0,1,  1,1,
		    -1,-1, 0,  0,0,1,  0,0,
		     1,-1, 0,  0,0,1,  1,0 ];

   S.octa = [      
            1,0,0,   1,1,1,   0,0,
            0,1,0,   1,1,1,   1,0,
            0,0,1,   1,1,1,   0,1,
      
            1,0,0,   1,-1,1,   0,0,
            0,-1,0,   1,-1,1,   1,0,
            0,0,1,   1,-1,1,   0,1,
      
            -1,0,0,   -1,1,1,   0,0,
            0,1,0,   -1,1,1,   1,0,
            0,0,1,   -1,1,1,   0,1,
      
            -1,0,0,   -1,-1,1,   0,0,
            0,-1,0,   -1,-1,1,   1,0,
            0,0,1,   -1,-1,1,   0,1,
      
            1,0,0,   1,1,-1,   0,0,
            0,1,0,   1,1,-1,   1,0,
            0,0,-1,   1,1,-1,   0,1,
      
            1,0,0,   1,-1,-1,   0,0,
            0,-1,0,   1,-1,-1,   1,0,
            0,0,-1,   1,-1,-1,   0,1,
      
            -1,0,0,   -1,1,-1,   0,0,
            0,1,0,   -1,1,-1,   1,0,
            0,0,-1,   -1,1,-1,   0,1,
      
            -1,0,0,   -1,-1,-1,   0,0,
            0,-1,0,   -1,-1,-1,   1,0,
            0,0,-1,   -1,-1,-1,   0,1,
      ];

   // GLUE TOGETHER TWO MESHES TO CREATE A SINGLE MESH

   let glueMeshes = (a,b) => {
      let mesh = a.slice();
      mesh.push(a.slice(a.length - S.VERTEX_SIZE, a.length));
      mesh.push(b.slice(0, S.VERTEX_SIZE));
      mesh.push(b);
      return mesh.flat();
   }

   // GIVEN A FUNCTION THAT MAPS (u,v) TO point AND normal,
   // AND GIVEN A MESH RESOLUTION, CREATE A PARAMETRIC MESH

   let uvMesh = (f,nu,nv) => {
      let mesh = [];
      for (let iv = 0 ; iv < nv ; iv++) {
         let v = iv / nv;
	 let strip = [];
         for (let iu = 0 ; iu <= nu ; iu++) {
	    let u = iu / nu;
	    strip = strip.concat(f(u,v));
	    strip = strip.concat(f(u,v+1/nv));
	 }
	 mesh = glueMeshes(mesh, strip);
      }
      return mesh;
   }

   // CREATE A UNIT SPHERE PARAMETRIC MESH

   S.coneMesh = uvMesh( (u,v) => {
      let theta = 2 * Math.PI * u;
      let rd = v;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      return [rd * cu, rd * su, rd,
              cu , su, 1,
              u, v];
 }, 20, 2);
  
  S.circleCap = uvMesh( (u,v) => {
      let theta = 2 * Math.PI * u;
      let rd = v;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      return [rd * cu, rd * su, 1,
              0, 0, 1,
              u, v];
  }, 20, 2) ;

  S.coneMesh = glueMeshes(S.coneMesh, S.circleCap );


   S.sphereMesh = uvMesh((u,v) => {
      let theta = 2 * Math.PI * u;
      let phi = Math.PI * v - Math.PI/2;
      let cu = Math.cos(theta);
      let su = Math.sin(theta);
      let cv = Math.cos(phi);
      let sv = Math.sin(phi);
      return [cu * cv, su * cv, sv,
              cu * cv, su * cv, sv,
	      u, v];
   }, 20, 10);

   // TRANSFORM A MESH BY A MATRIX ON THE CPU

   let transformMesh = (mesh, matrix) => {
      let result = [];
      let IMT = matrixTranspose(matrixInverse(matrix));
      for (let n = 0 ; n < mesh.length ; n += S.VERTEX_SIZE) {
         let V = mesh.slice(n, n + S.VERTEX_SIZE);
	 let P  = V.slice(0, 3);
	 let N  = V.slice(3, 6);
	 let UV = V.slice(6, 8);
	 P = matrixTransform(matrix, [P[0], P[1], P[2], 1]);
	 N = matrixTransform(IMT,    [N[0], N[1], N[2], 0]);
         result.push(P[0],P[1],P[2], N[0],N[1],N[2], UV);
      }
      return result.flat();
   }

   // A CUBE MESH IS SIX TRANSFORMED SQUARE MESHES GLUED TOGETHER

   let face0 = transformMesh(S.squareMesh, matrixTranslate([0,0,1]));
   let face1 = transformMesh(face0,        matrixRotx( Math.PI/2));
   let face2 = transformMesh(face0,        matrixRotx( Math.PI  ));
   let face3 = transformMesh(face0,        matrixRotx(-Math.PI/2));
   let face4 = transformMesh(face0,        matrixRoty(-Math.PI/2));
   let face5 = transformMesh(face0,        matrixRoty( Math.PI/2));
   S.cubeMesh = glueMeshes(face0,
                glueMeshes(face1,
                glueMeshes(face2,
                glueMeshes(face3,
                glueMeshes(face4,
		           face5)))));

    let octa0 = transformMesh(S.octa, matrixTranslate([0,0,0]));
    let octa1 = transformMesh(octa0,        matrixRotx( Math.PI/2));
    let octa2 = transformMesh(octa0,        matrixRotx( Math.PI  ));
    let octa3 = transformMesh(octa0,        matrixRotx(-Math.PI/2));
    let octa4 = transformMesh(octa0,        matrixRoty(-Math.PI/2));
    let octa5 = transformMesh(octa0,        matrixRoty( Math.PI/2));

    S.octaMesh = octa0.concat(octa1.concat(octa2.concat(octa3.concat(octa4.concat(octa5)))));

   // DRAW A SINGLE MESH. WE STILL NEED TO ADD MATERIAL PROPERTIES!

   S.drawMesh = (mesh, matrix) => {
      let gl = S.gl;
      S.setUniform('Matrix4fv', 'uMatrix', false, matrix);
      S.setUniform('Matrix4fv', 'uInvMatrix', false, matrixInverse(matrix));
      S.gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh), gl.STATIC_DRAW);
      S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, mesh.length / S.VERTEX_SIZE);
   }
   S.drawMeshcolor = (mesh, matrix, type, phong=S.whitePlastic) => {
      let gl = S.gl;
      S.setUniform('Matrix4fv', 'uPhong', false, phong.flat());
      S.setUniform('Matrix4fv', 'uMatrix', false, matrix);
      S.setUniform('Matrix4fv', 'uInvMatrix', false, matrixInverse(matrix));
      S.gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(mesh), gl.STATIC_DRAW);
      S.gl.drawArrays(S.gl.TRIANGLE_STRIP, 0, mesh.length / S.VERTEX_SIZE);
   }

`,
fragment: `
S.setFragmentShader(\`
   varying vec3 vPos, vNor;
   uniform float uOpacity;
   const int nL = \` + S.nL + \`;
   uniform vec3 uLd[nL];
   uniform vec3 uLc[nL];
   uniform mat4 uPhong;
   
vec3 shadeSurface(vec3 P, vec3 N, mat4 M){

   // EXTRACT PHONG PARAMETERS FROM MATERIAL MATRIX

   vec3  ambient  = M[0].rgb;
   vec3  diffuse  = M[1].rgb;
   vec3  specular = M[2].rgb;
   float p        = M[2].a;

   // COMPUTE NORMAL, INIT COLOR, APPROXIMATE VECTOR TO EYE

   vec3 c = ambient;
   vec3 E = vec3(0.,0.,1.);

   // LOOP THROUGH LIGHT SOURCES

   for (int l = 0 ; l < nL ; l++) {

     
      
         
      // COMPUTE DIFFUSE AND SPECULAR FOR THIS LIGHT SOURCE(reflection of light)

      
     vec3 R = 2. * dot(N, uLd[l]) * N - uLd[l];
     c += uLc[l] * (diffuse * max(0.,dot(N, uLd[l]))
                     + specular * pow(max(0., dot(R, E)), p));
   }
   //c *= 1. + .5 * noise(3.*N); // OPTIONAL SPOTTY TEXTURE
   return c;    
}

void main() {
      
   //float c = .2 + .8 * max(0.,dot(vNor,vec3(.57)));
   vec3 c = shadeSurface(vPos, vNor, uPhong);
   gl_FragColor = vec4(c,1.);
}
\`);
`,
vertex: `
S.setVertexShader(\`

   attribute vec3 aPos, aNor;
   varying   vec3 vPos, vNor;
   uniform   mat4 uMatrix, uInvMatrix, uProject;

   void main() {
      vec4 pos = uProject * uMatrix * vec4(aPos, 1.);
      vec4 nor = vec4(aNor, 0.) * uInvMatrix;
      vPos = pos.xyz;
      vNor = normalize(nor.xyz);
      gl_Position = pos * vec4(1.,1.,-.01,1.);
   }
\`)
`,
render: `

   // SET THE PROJECTION MATRIX BASED ON CAMERA FOCAL LENGTH

   let fl = 5.0;
   S.setUniform('Matrix4fv', 'uProject', false,
      [1,0,0,0, 0,1,0,0, 0,0,1,-1/fl, 0,0,0,1]);
   
   let add = (a,b) => [ a[0]+b[0], a[1]+b[1], a[2]+b[2] ];
   let dot = (a,b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
   let norm = v => Math.sqrt(dot(v,v));
   let normalize = v => { let s = norm(v); return [ v[0]/s, v[1]/s, v[2]/s ]; }
   let scale = (v,s) => [ s * v[0], s * v[1], s * v[2] ];
   let subtract = (a,b) => [ a[0]-b[0], a[1]-b[1], a[2]-b[2] ];
   let ldData = [ normalize([1,1,1]),
   normalize([-1,-1,-1]) ];
S.setUniform('3fv', 'uLd', ldData.flat());
S.setUniform('3fv', 'uLc', [ 1,1,1, .5,.3,.1 ]);

// DEFINE NUMBER OF LIGHTS FOR GPU

S.nL = ldData.length;
phong = S.redPlastic;

S.setUniform('Matrix4fv', 'uPhong', false, phong.flat());



   let m = new Matrix();

   // GET VALUES FROM THE HTML SLIDERS

   let T = 2 * time * rate.value / 100;
   let AL = .1 + .9 * arm_length.value / 100;
   let LL = .1 + .9 * leg_length.value / 100;
   let bodyHeight = .6;
   let pelvisHeight = .15;
   let neckHeight = .1;
   let shoulderLength = .2;

   // RENDER THE SCENE

   m.identity();
   m.roty(Math.sin(.5 * T));
   m.save();


      //PELVIS
      m.save();
        
        m.scale(.14,pelvisHeight/2,.07);
        S.drawMeshcolor(S.sphereMesh, m.get(),1,S.bluePlastic);
      m.restore();

      // LEGS

      for (let s = -1 ; s <= 1 ; s += 2) {
         let t = T - s * Math.PI/2;
         m.save();
           m.translate(s*.1,-.01 ,0);
           m.rotx(Math.cos(t));
           m.save();
              m.translate(0,-LL/2,0);
              m.scale(.05,LL/2,.05);
              S.drawMeshcolor(S.sphereMesh, m.get(),1,S.redPlastic);
           m.restore();
           m.translate(0,-LL,0);
           m.rotx(1 + Math.sin(t));
           m.save();
              m.translate(0,-LL/2,0);
              m.scale(.05,LL/2,.05);
              S.drawMeshcolor(S.sphereMesh, m.get(),1,S.redPlastic);
           m.restore();
           m.translate(0,-LL,0);
           m.rotx(-1.6-Math.sin(t));
           m.save();
              m.translate(0,-0.09,0);
              m.scale(.05,0.09,.05);
              S.drawMeshcolor(S.sphereMesh, m.get(),1,S.whitePlastic);
           m.restore();
         m.restore();
      }


      // CHEST
      m.save();
      m.translate(0,pelvisHeight/2, 0);
      m.roty(Math.cos(2.1*T)/4);
      m.save();
            
            m.translate(0,bodyHeight/2,0);
            m.scale(0.15, bodyHeight/2, 0.09);
            S.drawMeshcolor(S.sphereMesh, m.get(),1,S.orangePlastic);
        
         m.restore();
         // NECK
         m.translate(0,bodyHeight, 0);
         m.save();
               m.translate(0,neckHeight/2,0);
               m.scale(.025,neckHeight/2,.025);     
               S.drawMeshcolor(S.sphereMesh, m.get(),1,S.whitePlastic);
         m.restore();
         // HEAD
            m.save();
                  m.translate(0,neckHeight,0);
                  m.rotz(Math.cos(T));
                  m.save();
                     m.translate(0,.1,0);
                     m.scale(.1,.1,.1);
                     S.drawMeshcolor(S.sphereMesh, m.get(),1,S.purplePlastic);
                  m.restore();
               m.restore();
         //SHOULDER
         for (let s = -1 ; s <= 1 ; s += 2){
            let t = T + s * Math.PI/2;
            m.save();
                m.rotx(t);
                m.rotx(-s*Math.PI/12);
                    m.save();
                    m.translate(s*shoulderLength/2,0,0);
                    
                        m.scale(shoulderLength/2,.025,.025);
                        S.drawMeshcolor(S.sphereMesh, m.get(),1,S.redPlastic);
                    m.restore();
         //AMRS
            m.save()
            m.translate(s*shoulderLength,0,0);
            m.rotx(Math.cos(t));
            m.rotz(s);
            m.save();
              m.translate(0,-AL/2,0);
              m.scale(.035,AL/2,.035);
              S.drawMeshcolor(S.sphereMesh, m.get(),1,S.redPlastic);
           m.restore();
           m.translate(0,-AL,0);
           m.rotx(-1. + Math.sin(t));
           m.save();
              m.translate(0,-AL/2,0);
              m.scale(.035,AL/2,.035);
              S.drawMeshcolor(S.sphereMesh, m.get(),1,S.redPlastic);
              m.restore();
            //HAND
            m.save();
            m.translate(0,-AL,0);
            m.scale(.05,0.05,.05);
            S.drawMeshcolor(S.sphereMesh, m.get(),1,S.whitePlastic);
            m.restore();
            

            //WEAPON
            m.save();
               m.translate(0,-AL-0.025,0);
               m.scale(.2,.1,.2);
               m.rotx(-Math.PI/2)
               m.roty(Math.PI);
               S.drawMesh(S.coneMesh, m.get(), 0);

            m.restore();
            m.save();
   //OCTAHEDRON 
   m.translate(0,-AL-0.25,0);
   m.rotx(1.3*T);
   m.roty(2*T);
   m.save();
       m.scale(.1,.1,.1);
       S.drawMesh(S.octaMesh, m.get(), 1, S.whitePlastic);
   m.restore();
           m.restore();
         m.restore();



         }
      
   m.restore();

   
   
`,
events: `
   ;
`
};

}

