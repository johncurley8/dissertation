/*
 * John Curley, April 2025
 */

@#define N = 250
var 
    % Aggregate Variables
    c ${c}$ (long_name='Consumption')
    n ${n}$ (long_name='Total Labor')
    w ${w}$ (long_name='Wage')
    k ${k}$ (long_name='Total Capital')
    R ${R}$ (long_name='Effective Gross Interest Rate')
    d ${d}$ (long_name='Total Dividend')
    b ${b}$ (long_name='Total Bonds')

    % Firm-Specific Variables
@# for i in 1:N
    n@{i} ${n_@{i}}$ (long_name='Labor Type @{i}')
    k@{i} ${k_@{i}}$ (long_name='Capital Type @{i}')
    b@{i} ${b_@{i}}$ (long_name='Debt Type @{i}')
    d@{i} ${d_@{i}}$ (long_name='Dividend Type @{i}')
    y@{i} ${y_@{i}}$ (long_name='Output Type @{i}')
    mu@{i} ${\mu_@{i}}$ (long_name='Lagrange Multiplier Type @{i}')
    xi@{i} ${z_@{i}}$ (long_name='Financial Conditions Type @{i}')

    yhat@{i} ${\hat y@{i}}$ (long_name='Output deviation from trend')
    byhat@{i} ${\hat y}$ (long_name='Debt-repurchase to GDP ratio deviation from trend') 
    dyhat@{i} ${\hat y}$ (long_name='Dividend to GDP ratio deviation from trend')
    invest@{i} ${i@{i}}$ (long_name='Investment Type @{i}')
    ihat@{i} ${\hat i@{i}}$ (long_name='Investment deviation from trend')
    nhat@{i} ${\hat n@{i}}$ (long_name='Hours deviation from trend')
@# endfor

    % Quasi-exogenous (aggregate) variables
    z ${z}$ (long_name='Aggregate Technology ')
    xi ${\xi}$ (long_name='Financial Conditions')
    
    % Descriptive Variables
    y ${y}$ (long_name='Output')
    yhat ${\hat y}$ (long_name='output deviation from trend')
    chat ${\hat c}$ (long_name='consumption deviation from trend') 
    byhat ${\hat y}$ (long_name='debt-repurchase to GDP ratio deviation from trend') 
    dyhat ${\hat y}$ (long_name='dividend to GDP ratio deviation from trend')
    v ${v}$ (long_name='Value of the Firm')
    vyhat ${\hat y}$ (long_name='Equity value deviation from trend')
    r ${r}$ (long_name='Gross Interest Rate')
    invest ${i}$ (long_name='Investment')
    ihat ${\hat i}$ (long_name='investment deviation from trend')
    nhat ${\hat n}$ (long_name='hours deviation from trend')
    % muhat ${\hat \mu}$ (long_name='Lagrange multiplier from trend')
        ;

varexo eps_z ${\varepsilon_z}$ (long_name='Aggregate technology shock')
       eps_xi ${\varepsilon_{\xi}}$ (long_name='Aggregate financial Shock')
@# for i in 1:N
    eps_xi@{i} ${\varepsilon_{xi@{i}}}$ (long_name='Financial conditions type @{i} shock')
@# endfor
         ;

parameters theta ${\theta}$ (long_name='capital share')
        nu ${\nu}$ (long_name='labor share')
        betta ${\beta}$ (long_name='discount factor')
        alppha ${\alpha}$ (long_name='disutility from work')
        delta ${\delta}$ (long_name='depreciation')
        tau ${\tau}$ (long_name='tax wedge')
        kappa ${\kappa}$ (long_name='equity cost')
        siggma ${\sigma}$ (long_name='risk aversion')
        sigma_z ${\sigma_z}$ (long_name='std_z')
        sigma_xi ${\sigma_xi}$ (long_name='std_xi')
        @# for i in 1:N
            sigma_xi@{i}  % Firm-specific financial shock volatility
        @# endfor
        BY_ratio ${(\bar b/(1+\bar r)/\bar Y}$ (long_name='Debt output ratio')
        ppsi ${\psi}$ (long_name='psi')
        @# for i in 1:N
            rho_xi@{i}
        @# endfor
        % Aggregate shocks VAR matrix
        A11 ${A_{11}}$ (long_name='A_11')
        A12 ${A_{12}}$ (long_name='A_12')
        A21 ${A_{21}}$ (long_name='A_21')
        A22 ${A_{22}}$ (long_name='A_22')
        ;
     
    siggma=1;
    theta = 0.25;
    nu = 0.6;
    betta = 0.9825;
    delta = 0.025;
    tau = 0.35;        
    BY_ratio=3.36;
    ppsi = 1.0;

    kappa = 0.146;
    sigma_xi = 0.0098;
    @# for i in 1:N
        sigma_xi@{i} = sigma_xi*10;  % Firm-specific financial shocks
    @# endfor
    sigma_z = 0.0045;
    @# for i in 1:N
        rho_xi@{i} = 0.9;
    @# endfor
    A11 = 0.9457;
    A12 = 0;
    A21 = 0;
    A22 = 0.09703;
    options_.TeX=1;

model;
[name='FOC labor']
w/c^siggma - alppha/(1-n) = 0;

[name='Euler equation household']
c^(-siggma) = betta * ((R - tau) / (1 - tau)) * c(+1)^(-siggma);

[name='Budget constraint household']
w * n + b(-1) - b/R + d - c = 0;

@# for i in 1:N
    [name='FOC labor input firm @{i}']
    (nu) * z * k@{i}(-1)^theta * n@{i}^(nu-1) = w * (1 / (1 - mu@{i} * (1 + 2 * kappa * (d@{i} - steady_state(d@{i})))));
@# endfor

@# for i in 1:N
    [name='FOC capital firm @{i}']
    betta * (c / c(+1))^siggma * ((1 + 2*kappa*(d@{i} - steady_state(d@{i}))) / (1 + 2*kappa*(d@{i}(+1) - steady_state(d@{i})))) *
    (1 - delta + (1 - mu@{i}(+1) * (1 + 2*kappa*(d@{i}(+1) - steady_state(d@{i})))) * theta * z(+1) * k@{i}^(theta-1) * n@{i}(+1)^(nu))
    + xi * xi@{i} * mu@{i} * (1 + 2*kappa*(d@{i} - steady_state(d@{i}))) = 1;
@# endfor

@# for i in 1:N
    [name='FOC bonds firm @{i}']
    R * betta * (c / c(+1))^siggma * ((1 + 2*kappa*(d@{i} - steady_state(d@{i}))) / (1 + 2*kappa*(d@{i}(+1) - steady_state(d@{i})))) 
    + xi * xi@{i} * mu@{i} * (1 + 2*kappa*(d@{i} - steady_state(d@{i}))) * (R * (1 - tau) / (R - tau)) = 1;
@# endfor

@# for i in 1:N
    [name='Enforcement constraint firm @{i}']
    xi * xi@{i} *(k@{i} - b@{i} * ((1 - tau) / (R - tau))) = ppsi * z * k@{i}(-1)^theta * n@{i}^nu;
@# endfor

/* alternative enforcement constraint (section 6.1)
@# for i in 1:N
    [name='FOC bonds firm @{i}']
    R * betta * (c / c(+1))^siggma * ((1 + 2*kappa*(d@{i} - steady_state(d@{i}))) / (1 + 2*kappa*(d@{i}(+1) - steady_state(d@{i})))) 
    + mu@{i} * (1 + 2*kappa*(d@{i} - steady_state(d@{i}))) * (R * (1 - tau) / (R - tau)) = 1;
@# endfor

@# for i in 1:N
    [name='Enforcement constraint firm @{i}']
    xi * xi@{i} * k@{i} = b@{i} * ((1 - tau) / (R - tau)) + ppsi * z * k@{i}(-1)^theta * n@{i}^nu;
@# endfor
*/

@# for i in 1:N
    [name='Budget constraint firm @{i}']
    (1 - delta) * k@{i}(-1) + z * k@{i}(-1)^theta * n@{i}^nu - w * n@{i} - b@{i}(-1) + b@{i} / R - k@{i} - (d@{i} + kappa * (d@{i} - steady_state(d@{i}))^2) = 0;
@# endfor

@# for i in 1:N
    [name='Idiosyncratic financial conditions firm @{i}']
    log(xi@{i} / steady_state(xi@{i})) = rho_xi@{i} * log(xi@{i}(-1) / steady_state(xi@{i})) + eps_xi@{i};
@# endfor

@# for i in 1:N
    [name='Production function firm @{i}']
    y@{i} = z * k@{i}(-1)^theta * n@{i}^nu;
@# endfor

[name='Aggregate TFP process']
log(z/steady_state(z))=A11*log(z(-1)/steady_state(z))+A12*log(xi(-1)/steady_state(xi))+eps_z;

[name='Aggregate financial shock']
log(xi/steady_state(xi))=A22*log(xi(-1)/steady_state(xi))+eps_xi;


%%%% MARKET CLEARING CONDITIONS (APPENDED BY mkt_clearing.m) %%%%

[name='Capital market clearing']
k = k1 + k2 + k3 + k4 + k5 + k6 + k7 + k8 + k9 + k10 + k11 + k12 + k13 + k14 + k15 + k16 + k17 + k18 + k19 + k20 + k21 + k22 + k23 + k24 + k25 + k26 + k27 + k28 + k29 + k30 + k31 + k32 + k33 + k34 + k35 + k36 + k37 + k38 + k39 + k40 + k41 + k42 + k43 + k44 + k45 + k46 + k47 + k48 + k49 + k50 + k51 + k52 + k53 + k54 + k55 + k56 + k57 + k58 + k59 + k60 + k61 + k62 + k63 + k64 + k65 + k66 + k67 + k68 + k69 + k70 + k71 + k72 + k73 + k74 + k75 + k76 + k77 + k78 + k79 + k80 + k81 + k82 + k83 + k84 + k85 + k86 + k87 + k88 + k89 + k90 + k91 + k92 + k93 + k94 + k95 + k96 + k97 + k98 + k99 + k100 + k101 + k102 + k103 + k104 + k105 + k106 + k107 + k108 + k109 + k110 + k111 + k112 + k113 + k114 + k115 + k116 + k117 + k118 + k119 + k120 + k121 + k122 + k123 + k124 + k125 + k126 + k127 + k128 + k129 + k130 + k131 + k132 + k133 + k134 + k135 + k136 + k137 + k138 + k139 + k140 + k141 + k142 + k143 + k144 + k145 + k146 + k147 + k148 + k149 + k150 + k151 + k152 + k153 + k154 + k155 + k156 + k157 + k158 + k159 + k160 + k161 + k162 + k163 + k164 + k165 + k166 + k167 + k168 + k169 + k170 + k171 + k172 + k173 + k174 + k175 + k176 + k177 + k178 + k179 + k180 + k181 + k182 + k183 + k184 + k185 + k186 + k187 + k188 + k189 + k190 + k191 + k192 + k193 + k194 + k195 + k196 + k197 + k198 + k199 + k200 + k201 + k202 + k203 + k204 + k205 + k206 + k207 + k208 + k209 + k210 + k211 + k212 + k213 + k214 + k215 + k216 + k217 + k218 + k219 + k220 + k221 + k222 + k223 + k224 + k225 + k226 + k227 + k228 + k229 + k230 + k231 + k232 + k233 + k234 + k235 + k236 + k237 + k238 + k239 + k240 + k241 + k242 + k243 + k244 + k245 + k246 + k247 + k248 + k249 + k250;

[name='Labor market clearing']
n = n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8 + n9 + n10 + n11 + n12 + n13 + n14 + n15 + n16 + n17 + n18 + n19 + n20 + n21 + n22 + n23 + n24 + n25 + n26 + n27 + n28 + n29 + n30 + n31 + n32 + n33 + n34 + n35 + n36 + n37 + n38 + n39 + n40 + n41 + n42 + n43 + n44 + n45 + n46 + n47 + n48 + n49 + n50 + n51 + n52 + n53 + n54 + n55 + n56 + n57 + n58 + n59 + n60 + n61 + n62 + n63 + n64 + n65 + n66 + n67 + n68 + n69 + n70 + n71 + n72 + n73 + n74 + n75 + n76 + n77 + n78 + n79 + n80 + n81 + n82 + n83 + n84 + n85 + n86 + n87 + n88 + n89 + n90 + n91 + n92 + n93 + n94 + n95 + n96 + n97 + n98 + n99 + n100 + n101 + n102 + n103 + n104 + n105 + n106 + n107 + n108 + n109 + n110 + n111 + n112 + n113 + n114 + n115 + n116 + n117 + n118 + n119 + n120 + n121 + n122 + n123 + n124 + n125 + n126 + n127 + n128 + n129 + n130 + n131 + n132 + n133 + n134 + n135 + n136 + n137 + n138 + n139 + n140 + n141 + n142 + n143 + n144 + n145 + n146 + n147 + n148 + n149 + n150 + n151 + n152 + n153 + n154 + n155 + n156 + n157 + n158 + n159 + n160 + n161 + n162 + n163 + n164 + n165 + n166 + n167 + n168 + n169 + n170 + n171 + n172 + n173 + n174 + n175 + n176 + n177 + n178 + n179 + n180 + n181 + n182 + n183 + n184 + n185 + n186 + n187 + n188 + n189 + n190 + n191 + n192 + n193 + n194 + n195 + n196 + n197 + n198 + n199 + n200 + n201 + n202 + n203 + n204 + n205 + n206 + n207 + n208 + n209 + n210 + n211 + n212 + n213 + n214 + n215 + n216 + n217 + n218 + n219 + n220 + n221 + n222 + n223 + n224 + n225 + n226 + n227 + n228 + n229 + n230 + n231 + n232 + n233 + n234 + n235 + n236 + n237 + n238 + n239 + n240 + n241 + n242 + n243 + n244 + n245 + n246 + n247 + n248 + n249 + n250;

[name='Bond market clearing']
b = b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10 + b11 + b12 + b13 + b14 + b15 + b16 + b17 + b18 + b19 + b20 + b21 + b22 + b23 + b24 + b25 + b26 + b27 + b28 + b29 + b30 + b31 + b32 + b33 + b34 + b35 + b36 + b37 + b38 + b39 + b40 + b41 + b42 + b43 + b44 + b45 + b46 + b47 + b48 + b49 + b50 + b51 + b52 + b53 + b54 + b55 + b56 + b57 + b58 + b59 + b60 + b61 + b62 + b63 + b64 + b65 + b66 + b67 + b68 + b69 + b70 + b71 + b72 + b73 + b74 + b75 + b76 + b77 + b78 + b79 + b80 + b81 + b82 + b83 + b84 + b85 + b86 + b87 + b88 + b89 + b90 + b91 + b92 + b93 + b94 + b95 + b96 + b97 + b98 + b99 + b100 + b101 + b102 + b103 + b104 + b105 + b106 + b107 + b108 + b109 + b110 + b111 + b112 + b113 + b114 + b115 + b116 + b117 + b118 + b119 + b120 + b121 + b122 + b123 + b124 + b125 + b126 + b127 + b128 + b129 + b130 + b131 + b132 + b133 + b134 + b135 + b136 + b137 + b138 + b139 + b140 + b141 + b142 + b143 + b144 + b145 + b146 + b147 + b148 + b149 + b150 + b151 + b152 + b153 + b154 + b155 + b156 + b157 + b158 + b159 + b160 + b161 + b162 + b163 + b164 + b165 + b166 + b167 + b168 + b169 + b170 + b171 + b172 + b173 + b174 + b175 + b176 + b177 + b178 + b179 + b180 + b181 + b182 + b183 + b184 + b185 + b186 + b187 + b188 + b189 + b190 + b191 + b192 + b193 + b194 + b195 + b196 + b197 + b198 + b199 + b200 + b201 + b202 + b203 + b204 + b205 + b206 + b207 + b208 + b209 + b210 + b211 + b212 + b213 + b214 + b215 + b216 + b217 + b218 + b219 + b220 + b221 + b222 + b223 + b224 + b225 + b226 + b227 + b228 + b229 + b230 + b231 + b232 + b233 + b234 + b235 + b236 + b237 + b238 + b239 + b240 + b241 + b242 + b243 + b244 + b245 + b246 + b247 + b248 + b249 + b250;

[name='Dividend market clearing']
d = d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9 + d10 + d11 + d12 + d13 + d14 + d15 + d16 + d17 + d18 + d19 + d20 + d21 + d22 + d23 + d24 + d25 + d26 + d27 + d28 + d29 + d30 + d31 + d32 + d33 + d34 + d35 + d36 + d37 + d38 + d39 + d40 + d41 + d42 + d43 + d44 + d45 + d46 + d47 + d48 + d49 + d50 + d51 + d52 + d53 + d54 + d55 + d56 + d57 + d58 + d59 + d60 + d61 + d62 + d63 + d64 + d65 + d66 + d67 + d68 + d69 + d70 + d71 + d72 + d73 + d74 + d75 + d76 + d77 + d78 + d79 + d80 + d81 + d82 + d83 + d84 + d85 + d86 + d87 + d88 + d89 + d90 + d91 + d92 + d93 + d94 + d95 + d96 + d97 + d98 + d99 + d100 + d101 + d102 + d103 + d104 + d105 + d106 + d107 + d108 + d109 + d110 + d111 + d112 + d113 + d114 + d115 + d116 + d117 + d118 + d119 + d120 + d121 + d122 + d123 + d124 + d125 + d126 + d127 + d128 + d129 + d130 + d131 + d132 + d133 + d134 + d135 + d136 + d137 + d138 + d139 + d140 + d141 + d142 + d143 + d144 + d145 + d146 + d147 + d148 + d149 + d150 + d151 + d152 + d153 + d154 + d155 + d156 + d157 + d158 + d159 + d160 + d161 + d162 + d163 + d164 + d165 + d166 + d167 + d168 + d169 + d170 + d171 + d172 + d173 + d174 + d175 + d176 + d177 + d178 + d179 + d180 + d181 + d182 + d183 + d184 + d185 + d186 + d187 + d188 + d189 + d190 + d191 + d192 + d193 + d194 + d195 + d196 + d197 + d198 + d199 + d200 + d201 + d202 + d203 + d204 + d205 + d206 + d207 + d208 + d209 + d210 + d211 + d212 + d213 + d214 + d215 + d216 + d217 + d218 + d219 + d220 + d221 + d222 + d223 + d224 + d225 + d226 + d227 + d228 + d229 + d230 + d231 + d232 + d233 + d234 + d235 + d236 + d237 + d238 + d239 + d240 + d241 + d242 + d243 + d244 + d245 + d246 + d247 + d248 + d249 + d250;

[name='Output market clearing']
y = y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8 + y9 + y10 + y11 + y12 + y13 + y14 + y15 + y16 + y17 + y18 + y19 + y20 + y21 + y22 + y23 + y24 + y25 + y26 + y27 + y28 + y29 + y30 + y31 + y32 + y33 + y34 + y35 + y36 + y37 + y38 + y39 + y40 + y41 + y42 + y43 + y44 + y45 + y46 + y47 + y48 + y49 + y50 + y51 + y52 + y53 + y54 + y55 + y56 + y57 + y58 + y59 + y60 + y61 + y62 + y63 + y64 + y65 + y66 + y67 + y68 + y69 + y70 + y71 + y72 + y73 + y74 + y75 + y76 + y77 + y78 + y79 + y80 + y81 + y82 + y83 + y84 + y85 + y86 + y87 + y88 + y89 + y90 + y91 + y92 + y93 + y94 + y95 + y96 + y97 + y98 + y99 + y100 + y101 + y102 + y103 + y104 + y105 + y106 + y107 + y108 + y109 + y110 + y111 + y112 + y113 + y114 + y115 + y116 + y117 + y118 + y119 + y120 + y121 + y122 + y123 + y124 + y125 + y126 + y127 + y128 + y129 + y130 + y131 + y132 + y133 + y134 + y135 + y136 + y137 + y138 + y139 + y140 + y141 + y142 + y143 + y144 + y145 + y146 + y147 + y148 + y149 + y150 + y151 + y152 + y153 + y154 + y155 + y156 + y157 + y158 + y159 + y160 + y161 + y162 + y163 + y164 + y165 + y166 + y167 + y168 + y169 + y170 + y171 + y172 + y173 + y174 + y175 + y176 + y177 + y178 + y179 + y180 + y181 + y182 + y183 + y184 + y185 + y186 + y187 + y188 + y189 + y190 + y191 + y192 + y193 + y194 + y195 + y196 + y197 + y198 + y199 + y200 + y201 + y202 + y203 + y204 + y205 + y206 + y207 + y208 + y209 + y210 + y211 + y212 + y213 + y214 + y215 + y216 + y217 + y218 + y219 + y220 + y221 + y222 + y223 + y224 + y225 + y226 + y227 + y228 + y229 + y230 + y231 + y232 + y233 + y234 + y235 + y236 + y237 + y238 + y239 + y240 + y241 + y242 + y243 + y244 + y245 + y246 + y247 + y248 + y249 + y250;



/*
end market clearing
*/

[name='Law of Motion of Capital']
invest=k-(1-delta)*k(-1);

[name='Definition output deviations from trend']
yhat=log(y)-log(steady_state(y));

[name='Definition consumption deviation from trend']
chat=log(c)-log(steady_state(c));

[name='Definition debt repurchase share in GDP']
byhat=(b(-1)/(1+r(-1))-b/(1+r))/y;    

[name='Definition equity payout to GDP ratio']
dyhat=d/y;

[name='Definition firm value']
v=d+betta*c/c(+1)*v(+1);

[name='Definition equity share']
vyhat=log(v/(k(-1)-b(-1)))-log(steady_state(v/(k(-1)-b(-1))));

[name='Definition before tax interest rate']
r=(R-tau)/(1-tau)-1;

[name='Definition investment deviation from trend']
ihat=log(invest)-log(steady_state(invest));

[name='Definition hours deviation from trend']
nhat=log(n)-log(steady_state(n));

@# for i in 1:N
    [name='Firm 1 law of Motion of Capital']
    invest@{i}=k@{i}-(1-delta)*k@{i}(-1);
@# endfor

@# for i in 1:N
    [name='Firm 1 definition output deviations from trend']
    yhat@{i}=log(y@{i})-log(steady_state(y@{i}));
@# endfor

@# for i in 1:N
    [name='Firm 1 definition debt repurchase share in GDP']
    byhat@{i}=(b@{i}(-1)/(1+r(-1))-b@{i}/(1+r))/y@{i};   
@# endfor

@# for i in 1:N
    [name='Firm 1 equity payout to GDP ratio']
    dyhat@{i}=d@{i}/y@{i};
@# endfor

@# for i in 1:N
    [name='Definition investment deviation from trend']
    ihat@{i}=log(invest@{i})-log(steady_state(invest@{i}));
@# endfor

@# for i in 1:N
    [name='Firm 1 definition hours deviation from trend']
    nhat@{i}=log(n@{i})-log(steady_state(n@{i}));
@# endfor

end;


shocks;
    var eps_z= sigma_z^2;
    var eps_xi= sigma_xi^2;
@# for i in 1:N
    var eps_xi@{i} = sigma_xi@{i}^2;
@# endfor 
end;

options_.dynatol.f = 1e-2;  % Set function tolerance
options_.dynatol.x = 1e-2; 
steady;

write_latex_original_model;
write_latex_dynamic_model;
write_latex_parameter_table;

% make VAR diagonal for generation of IRFS
set_param_value('A12',0);
set_param_value('A21',0);


% Run first-order perturbation
stoch_simul(order=1, irf=50, nograph);
irfs_order1 = oo_.irfs; 
save('IRFs_het_xi.mat', 'oo_');

var_decomp = oo_.variance_decomposition;
shock_names = cellstr(M_.exo_names);
var_names   = cellstr(M_.endo_names);

save('var_decomp_xi.mat', 'var_decomp', 'shock_names', 'var_names');

/*
%% Display some statistics
fprintf('(b/(1+r)/Y: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))/(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact'))))/oo_.dr.ys(strmatch('y',M_.endo_names,'exact')))
fprintf('Debt-Capital Ratio: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))/(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact'))))/oo_.dr.ys(strmatch('k',M_.endo_names,'exact')))
fprintf('Total Debt-Output Ratio: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))+oo_.dr.ys(strmatch('y',M_.endo_names,'exact')))/(oo_.dr.ys(strmatch('y',M_.endo_names,'exact'))*(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact')))))
fprintf('Total Debt-Capital Ratio: %4.3f\n',(oo_.dr.ys(strmatch('b',M_.endo_names,'exact'))+oo_.dr.ys(strmatch('y',M_.endo_names,'exact')))/(oo_.dr.ys(strmatch('k',M_.endo_names,'exact'))*(1+oo_.dr.ys(strmatch('r',M_.endo_names,'exact')))))
fprintf('Capital-Output Ratio: %4.3f\n',(oo_.dr.ys(strmatch('k',M_.endo_names,'exact'))/oo_.dr.ys(strmatch('y',M_.endo_names,'exact'))))
*/

%%%%%%%%%%%%%%%%%%%%%%%% Create Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeRange = 0:30;

/*
figure(1);
set(gcf, 'Name', 'Figure 1: IRF Comparison (First vs. Second Order)', 'NumberTitle', 'off');

subplot(2,4,1)
plot(-irfs_order1.yhat_eps_z * 100, 'b', 'LineWidth', 1.5) % First-order
hold on
plot(-irfs_order2.yhat_eps_z * 100, 'r--', 'LineWidth', 1.5) % Second-order
title('Output')

subplot(2,4,2)
plot(-irfs_order1.chat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.chat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Consumption')

subplot(2,4,3)
plot(-irfs_order1.nhat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.nhat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Hours')

subplot(2,4,4)
plot(-irfs_order1.ihat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.ihat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Investment')

subplot(2,4,5)
plot(-irfs_order1.byhat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.byhat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Debt rep./Y')

subplot(2,4,6)
plot(-irfs_order1.dyhat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.dyhat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Equity Payout/Y')

subplot(2,4,7)
plot(-irfs_order1.vyhat_eps_z * 100, 'b', 'LineWidth', 1.5)
hold on
plot(-irfs_order2.vyhat_eps_z * 100, 'r--', 'LineWidth', 1.5)
title('Equity Value')

legend('First Order', 'Second Order')

figure(2);
set(gcf, 'Name', 'Figure 1: Aggregate Impulse Responses', 'NumberTitle', 'off');
        subplot(2,4,1)
        plot(-oo_.irfs.yhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.yhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -0.8 0]);
        ll = legend('TFP shock', 'Financial Shock', 'Location', 'southwest');
        title('Output')
        
        subplot(2,4,2)
        plot(-oo_.irfs.chat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.chat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -0.6 0.06]);
        title('Consumption')
        
        subplot(2,4,3)
        plot(-oo_.irfs.nhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.nhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -1.2 0.8]);
        title('Hours')
        
        subplot(2,4,4)
        plot(-oo_.irfs.ihat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.ihat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -3.5 0.5]);
        title('Investment')
  
        subplot(2,4,5)
        plot(-oo_.irfs.byhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.byhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -1.9 2.4]);
        title('Debt rep./Y')
        
        subplot(2,4,6)
        plot(-oo_.irfs.dyhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.dyhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -1.5 1.5]);
        title('Equity Payout/Y')
        
        subplot(2,4,7)
        plot(-oo_.irfs.vyhat_eps_z*100,'b')
        hold on
        plot(-oo_.irfs.vyhat_eps_xi*100,'r--')
        axis([min(timeRange) max(timeRange) -0.3 0.4]);
        title('Equity Value')
*/