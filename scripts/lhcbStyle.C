<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
  <meta http-equiv="content-type" content="text/html;charset=UTF-8" />
  <meta http-equiv="generator" content="WebSVN 2.3.3" /> <!-- leave this for stats -->
  <link rel="shortcut icon" type="image/x-icon" href="/websvn/templates/calm/images/favicon.ico" />
  <link type="text/css" href="/websvn/templates/calm/styles.css" rel="stylesheet" media="screen" />
  <link rel='alternate' type='application/rss+xml' title='WebSVN RSS' href='/world/wsvn/lhcb/Urania/trunk/RootTools/LHCbStyle/src/lhcbStyle.C?op=rss&amp;' />
  <!--[if gte IE 5.5000]>
  <script type="text/javascript" src="/websvn/templates/calm/png.js"></script>
  <style type="text/css" media="screen">
  tbody tr td { padding:1px 0 }
  #wrap h2 { padding:10px 5px 0 5px; margin-bottom:-8px }
  </style>
  <![endif]-->
  <title>
       WebSVN
          - lhcb
               - Rev 188463
            - /Urania/trunk/RootTools/LHCbStyle/src/lhcbStyle.C
  </title>
  <script type="text/javascript">
  //<![CDATA[
       function getPath()
       {
         return '/websvn';
       }
       
       function checkCB(chBox)
       {
          count = 0
          first = null
          f = chBox.form
          for (i = 0 ; i < f.elements.length ; i++)
          if (f.elements[i].type == 'checkbox' && f.elements[i].checked)
          {
             if (first == null && f.elements[i] != chBox)
                first = f.elements[i]
             count += 1
          }
          
          if (count > 2) 
          {
             first.checked = false
             count -= 1
          }
       }
  //]]>
  </script>
</head>
<body id="file">
<div id="container">
	<div id="select">
		<form method="get" action="" id="project"><input type="hidden" name="op" value="rep" /><select name="repname" onchange="javascript:this.form.submit();"><option value="lhcb" selected="selected">lhcb</option><option value="lhcb" selected="selected">lhcb</option></select><noscript><input type="submit" value="Go" /></noscript></form>
		<form method="get" action="" id="template"><select name="template" onchange="javascript:this.form.submit();"><option value="BlueGrey">BlueGrey</option><option value="calm" selected="selected">calm</option><option value="Elegant">Elegant</option></select><noscript><input type="submit" value="Go" /></noscript></form>
		<form method="get" action="" id="language"><select name="language" onchange="javascript:this.form.submit();"><option value="ca">Catal&agrave;-Valenci&agrave; - Catalan</option><option value="zh-CN">&#20013;&#25991; - Chinese (Simplified)</option><option value="zh-TW">&#20013;&#25991; - Chinese (Traditional)</option><option value="cs">&#268;esky - Czech</option><option value="da">Dansk - Danish</option><option value="nl">Nederlands - Dutch</option><option value="en" selected="selected">English - English</option><option value="fi">Suomi - Finnish</option><option value="fr">Fran&ccedil;ais - French</option><option value="de">Deutsch - German</option><option value="he-IL">&#1506;&#1489;&#1512;&#1497;&#1514; - Hebrew</option><option value="hin">&#2361;&#2367;&#2306;&#2342;&#2368; - Hindi</option><option value="hu">Magyar - Hungarian</option><option value="id">Bahasa Indonesia - Indonesian</option><option value="it">Italiano - Italian</option><option value="ja">&#26085;&#26412;&#35486; - Japanese</option><option value="ko">&#54620;&#44397;&#50612; - Korean</option><option value="mk">&#1052;&#1072;&#1082;&#1077;&#1076;&#1086;&#1085;&#1089;&#1082;&#1080; - Macedonian</option><option value="mr">&#2350;&#2352;&#2366;&#2336;&#2368; - Marathi</option><option value="no">Norsk - Norwegian</option><option value="pl">Polski - Polish</option><option value="pt">Portugu&ecirc;s - Portuguese</option><option value="pt-BR">Portugu&ecirc;s - Portuguese (Brazil)</option><option value="ru">&#1056;&#1091;&#1089;&#1089;&#1082;&#1080;&#1081; - Russian</option><option value="sk">Sloven&#269;ina - Slovak</option><option value="sl">Sloven&#353;&#269;ina - Slovenian</option><option value="es">Espa&ntilde;ol - Spanish</option><option value="sv">Svenska - Swedish</option><option value="tr">T&uuml;rk&ccedil;e - Turkish</option><option value="uk">&#1059;&#1082;&#1088;&#1072;&#1111;&#1085;&#1089;&#1100;&#1082;&#1072; - Ukrainian</option><option value="uz">O&euml;zbekcha - Uzbek</option></select><noscript><input type="submit" value="Go" /></noscript></form>
	</div>
	<h1><a href="/world/wsvn/?" title="Subversion Repositories">Subversion Repositories</a>
		<span><a href="/world/wsvn/lhcb?">lhcb</a></span>
	</h1>
  <h2 id="pathlinks"><a href="/world/wsvn/lhcb/?" class="root"><span>(root)</span></a>/<a href="/world/wsvn/lhcb/Urania/?#acef0b56f9ef253b0c2680e3064fd8ce3">Urania</a>/<a href="/world/wsvn/lhcb/Urania/trunk/?#a1148f2bfe82b0824ab7e4b7750adb162">trunk</a>/<a href="/world/wsvn/lhcb/Urania/trunk/RootTools/?#a0d34be8d8e24fef38658a785c3472860">RootTools</a>/<a href="/world/wsvn/lhcb/Urania/trunk/RootTools/LHCbStyle/?#ad69e9e51e54de48f5f407fe7750524ad">LHCbStyle</a>/<a href="/world/wsvn/lhcb/Urania/trunk/RootTools/LHCbStyle/src/?#ae0d1341bb4e971c72a9f0e572f694c9f">src</a>/<span class="file">lhcbStyle.C</span> - Rev 188463</h2>
  <div id="revjump"><form method="get" action="" id="revision">Rev <input type="text" size="5" name="rev" placeholder="HEAD" /><span class="submit"><input type="submit" value="Go" /></span></form></div>
  <p>
    <span class="prev"><a href="/world/wsvn/lhcb/Urania/trunk/RootTools/LHCbStyle/src/lhcbStyle.C?rev=154683">Rev 154683</a></span> &#124;
    <span class="blame"><a href="/world/wsvn/lhcb/Urania/trunk/RootTools/LHCbStyle/src/lhcbStyle.C?op=blame&amp;rev=188463">Blame</a></span> &#124;
    <span class="diff"><a href="/world/wsvn/lhcb/Urania/trunk/RootTools/LHCbStyle/src/lhcbStyle.C?op=diff&amp;rev=188463">Compare with Previous</a></span> &#124;
    <span class="changes"><a href="/world/wsvn/lhcb/Urania/trunk/RootTools/LHCbStyle/src/lhcbStyle.C?op=revision&amp;rev=188463">Last modification</a></span> &#124;
    <span class="log"><a href="/world/wsvn/lhcb/Urania/trunk/RootTools/LHCbStyle/src/lhcbStyle.C?op=log&amp;rev=188463">View Log</a></span>
    &#124; <span class="compress"><a href="/world/wsvn/lhcb/Urania/trunk/RootTools/LHCbStyle/src/lhcbStyle.C?op=dl&amp;rev=188463">Download</a></span>
    &#124; <span class="feed"><a href="/world/wsvn/lhcb/Urania/trunk/RootTools/LHCbStyle/src/lhcbStyle.C?op=rss&amp;">RSS feed</a></span>
  </p>
  <div class="listing">
<pre><code>// all users - please change the name of this file to lhcbStyle.C
</code><code>// Commits to lhcbdocs svn of .C files are not allowed
</code><code>{
</code><code>
</code><code>  // define names for colours
</code><code>  Int_t black  = 1;
</code><code>  Int_t red    = 2;
</code><code>  Int_t green  = 3;
</code><code>  Int_t blue   = 4;
</code><code>  Int_t yellow = 5; 
</code><code>  Int_t magenta= 6;
</code><code>  Int_t cyan   = 7;
</code><code>  Int_t purple = 9;
</code><code>  
</code><code>
</code><code>////////////////////////////////////////////////////////////////////
</code><code>// PURPOSE:
</code><code>//
</code><code>// This macro defines a standard style for (black-and-white) 
</code><code>// &quot;publication quality&quot; LHCb ROOT plots. 
</code><code>//
</code><code>// USAGE:
</code><code>//
</code><code>// Include the lines
</code><code>//   gROOT-&gt;ProcessLine(&quot;.L lhcbstyle.C&quot;);
</code><code>//   lhcbStyle();
</code><code>// at the beginning of your root macro.
</code><code>//
</code><code>// Example usage is given in myPlot.C
</code><code>//
</code><code>// COMMENTS:
</code><code>//
</code><code>// Font:
</code><code>// 
</code><code>// The font is chosen to be 132, this is Times New Roman (like the text of
</code><code>//  your document) with precision 2.
</code><code>//
</code><code>// &quot;Landscape histograms&quot;:
</code><code>//
</code><code>// The style here is designed for more or less square plots.
</code><code>// For longer histograms, or canvas with many pads, adjustements are needed. 
</code><code>// For instance, for a canvas with 1x5 histograms:
</code><code>//  TCanvas* c1 = new TCanvas(&quot;c1&quot;, &quot;L0 muons&quot;, 600, 800);
</code><code>//  c1-&gt;Divide(1,5);
</code><code>//  Adaptions like the following will be needed:
</code><code>//  gStyle-&gt;SetTickLength(0.05,&quot;x&quot;);
</code><code>//  gStyle-&gt;SetTickLength(0.01,&quot;y&quot;);
</code><code>//  gStyle-&gt;SetLabelSize(0.15,&quot;x&quot;);
</code><code>//  gStyle-&gt;SetLabelSize(0.1,&quot;y&quot;);
</code><code>//  gStyle-&gt;SetStatW(0.15);
</code><code>//  gStyle-&gt;SetStatH(0.5);
</code><code>//
</code><code>// Authors: Thomas Schietinger, Andrew Powell, Chris Parkes, Niels Tuning
</code><code>// Maintained by Editorial board member (currently Niels)
</code><code>///////////////////////////////////////////////////////////////////
</code><code>
</code><code>  // Use times new roman, precision 2 
</code><code>  Int_t lhcbFont        = 132;  // Old LHCb style: 62;
</code><code>  // Line thickness
</code><code>  Double_t lhcbWidth    = 2.00; // Old LHCb style: 3.00;
</code><code>  // Text size
</code><code>  Double_t lhcbTSize    = 0.06; 
</code><code>  
</code><code>  // use plain black on white colors
</code><code>  gROOT-&gt;SetStyle(&quot;Plain&quot;); 
</code><code>  TStyle *lhcbStyle= new TStyle(&quot;lhcbStyle&quot;,&quot;LHCb plots style&quot;);
</code><code>  
</code><code>  //lhcbStyle-&gt;SetErrorX(0); //  don&apos;t suppress the error bar along X
</code><code>
</code><code>  lhcbStyle-&gt;SetFillColor(1);
</code><code>  lhcbStyle-&gt;SetFillStyle(1001);   // solid
</code><code>  lhcbStyle-&gt;SetFrameFillColor(0);
</code><code>  lhcbStyle-&gt;SetFrameBorderMode(0);
</code><code>  lhcbStyle-&gt;SetPadBorderMode(0);
</code><code>  lhcbStyle-&gt;SetPadColor(0);
</code><code>  lhcbStyle-&gt;SetCanvasBorderMode(0);
</code><code>  lhcbStyle-&gt;SetCanvasColor(0);
</code><code>  lhcbStyle-&gt;SetStatColor(0);
</code><code>  lhcbStyle-&gt;SetLegendBorderSize(0);
</code><code>  lhcbStyle-&gt;SetLegendFont(132);
</code><code>
</code><code>  // If you want the usual gradient palette (blue -&gt; red)
</code><code>  lhcbStyle-&gt;SetPalette(1);
</code><code>  // If you want colors that correspond to gray scale in black and white:
</code><code>  int colors[8] = {0,5,7,3,6,2,4,1};
</code><code>  lhcbStyle-&gt;SetPalette(8,colors);
</code><code>
</code><code>  // set the paper &amp; margin sizes
</code><code>  lhcbStyle-&gt;SetPaperSize(20,26);
</code><code>  lhcbStyle-&gt;SetPadTopMargin(0.05);
</code><code>  lhcbStyle-&gt;SetPadRightMargin(0.05); // increase for colz plots
</code><code>  lhcbStyle-&gt;SetPadBottomMargin(0.16);
</code><code>  lhcbStyle-&gt;SetPadLeftMargin(0.14);
</code><code>  
</code><code>  // use large fonts
</code><code>  lhcbStyle-&gt;SetTextFont(lhcbFont);
</code><code>  lhcbStyle-&gt;SetTextSize(lhcbTSize);
</code><code>  lhcbStyle-&gt;SetLabelFont(lhcbFont,&quot;x&quot;);
</code><code>  lhcbStyle-&gt;SetLabelFont(lhcbFont,&quot;y&quot;);
</code><code>  lhcbStyle-&gt;SetLabelFont(lhcbFont,&quot;z&quot;);
</code><code>  lhcbStyle-&gt;SetLabelSize(lhcbTSize,&quot;x&quot;);
</code><code>  lhcbStyle-&gt;SetLabelSize(lhcbTSize,&quot;y&quot;);
</code><code>  lhcbStyle-&gt;SetLabelSize(lhcbTSize,&quot;z&quot;);
</code><code>  lhcbStyle-&gt;SetTitleFont(lhcbFont);
</code><code>  lhcbStyle-&gt;SetTitleFont(lhcbFont,&quot;x&quot;);
</code><code>  lhcbStyle-&gt;SetTitleFont(lhcbFont,&quot;y&quot;);
</code><code>  lhcbStyle-&gt;SetTitleFont(lhcbFont,&quot;z&quot;);
</code><code>  lhcbStyle-&gt;SetTitleSize(1.2*lhcbTSize,&quot;x&quot;);
</code><code>  lhcbStyle-&gt;SetTitleSize(1.2*lhcbTSize,&quot;y&quot;);
</code><code>  lhcbStyle-&gt;SetTitleSize(1.2*lhcbTSize,&quot;z&quot;);
</code><code>
</code><code>  // use medium bold lines and thick markers
</code><code>  lhcbStyle-&gt;SetLineWidth(lhcbWidth);
</code><code>  lhcbStyle-&gt;SetFrameLineWidth(lhcbWidth);
</code><code>  lhcbStyle-&gt;SetHistLineWidth(lhcbWidth);
</code><code>  lhcbStyle-&gt;SetFuncWidth(lhcbWidth);
</code><code>  lhcbStyle-&gt;SetGridWidth(lhcbWidth);
</code><code>  lhcbStyle-&gt;SetLineStyleString(2,&quot;[12 12]&quot;); // postscript dashes
</code><code>  lhcbStyle-&gt;SetMarkerStyle(20);
</code><code>  lhcbStyle-&gt;SetMarkerSize(1.0);
</code><code>
</code><code>  // label offsets
</code><code>  lhcbStyle-&gt;SetLabelOffset(0.010,&quot;X&quot;);
</code><code>  lhcbStyle-&gt;SetLabelOffset(0.010,&quot;Y&quot;);
</code><code>
</code><code>  // by default, do not display histogram decorations:
</code><code>  lhcbStyle-&gt;SetOptStat(0);  
</code><code>  //lhcbStyle-&gt;SetOptStat(&quot;emr&quot;);  // show only nent -e , mean - m , rms -r
</code><code>  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
</code><code>  lhcbStyle-&gt;SetStatFormat(&quot;6.3g&quot;); // specified as c printf options
</code><code>  lhcbStyle-&gt;SetOptTitle(0);
</code><code>  lhcbStyle-&gt;SetOptFit(0);
</code><code>  //lhcbStyle-&gt;SetOptFit(1011); // order is probability, Chi2, errors, parameters
</code><code>  //titles
</code><code>  lhcbStyle-&gt;SetTitleOffset(0.95,&quot;X&quot;);
</code><code>  lhcbStyle-&gt;SetTitleOffset(0.95,&quot;Y&quot;);
</code><code>  lhcbStyle-&gt;SetTitleOffset(1.2,&quot;Z&quot;);
</code><code>  lhcbStyle-&gt;SetTitleFillColor(0);
</code><code>  lhcbStyle-&gt;SetTitleStyle(0);
</code><code>  lhcbStyle-&gt;SetTitleBorderSize(0);
</code><code>  lhcbStyle-&gt;SetTitleFont(lhcbFont,&quot;title&quot;);
</code><code>  lhcbStyle-&gt;SetTitleX(0.0);
</code><code>  lhcbStyle-&gt;SetTitleY(1.0); 
</code><code>  lhcbStyle-&gt;SetTitleW(1.0);
</code><code>  lhcbStyle-&gt;SetTitleH(0.05);
</code><code>  
</code><code>  // look of the statistics box:
</code><code>  lhcbStyle-&gt;SetStatBorderSize(0);
</code><code>  lhcbStyle-&gt;SetStatFont(lhcbFont);
</code><code>  lhcbStyle-&gt;SetStatFontSize(0.05);
</code><code>  lhcbStyle-&gt;SetStatX(0.9);
</code><code>  lhcbStyle-&gt;SetStatY(0.9);
</code><code>  lhcbStyle-&gt;SetStatW(0.25);
</code><code>  lhcbStyle-&gt;SetStatH(0.15);
</code><code>
</code><code>  // put tick marks on top and RHS of plots
</code><code>  lhcbStyle-&gt;SetPadTickX(1);
</code><code>  lhcbStyle-&gt;SetPadTickY(1);
</code><code>
</code><code>  // histogram divisions: only 5 in x to avoid label overlaps
</code><code>  lhcbStyle-&gt;SetNdivisions(505,&quot;x&quot;);
</code><code>  lhcbStyle-&gt;SetNdivisions(510,&quot;y&quot;);
</code><code>  
</code><code>  gROOT-&gt;SetStyle(&quot;lhcbStyle&quot;);
</code><code>  gROOT-&gt;ForceStyle();
</code><code>
</code><code>  // add LHCb label
</code><code>  TPaveText* lhcbName = new TPaveText(gStyle-&gt;GetPadLeftMargin() + 0.05,
</code><code>                                      0.87 - gStyle-&gt;GetPadTopMargin(),
</code><code>                                      gStyle-&gt;GetPadLeftMargin() + 0.20,
</code><code>                                      0.95 - gStyle-&gt;GetPadTopMargin(),
</code><code>                                      &quot;BRNDC&quot;);
</code><code>  lhcbName-&gt;AddText(&quot;LHCb&quot;);
</code><code>  lhcbName-&gt;SetFillColor(0);
</code><code>  lhcbName-&gt;SetTextAlign(12);
</code><code>  lhcbName-&gt;SetBorderSize(0);
</code><code>
</code><code>  TText *lhcbLabel = new TText();
</code><code>  lhcbLabel-&gt;SetTextFont(lhcbFont);
</code><code>  lhcbLabel-&gt;SetTextColor(1);
</code><code>  lhcbLabel-&gt;SetTextSize(lhcbTSize);
</code><code>  lhcbLabel-&gt;SetTextAlign(12);
</code><code>
</code><code>  TLatex *lhcbLatex = new TLatex();
</code><code>  lhcbLatex-&gt;SetTextFont(lhcbFont);
</code><code>  lhcbLatex-&gt;SetTextColor(1);
</code><code>  lhcbLatex-&gt;SetTextSize(lhcbTSize);
</code><code>  lhcbLatex-&gt;SetTextAlign(12);
</code><code>
</code><code>  cout &lt;&lt; &quot;-------------------------&quot; &lt;&lt; endl;  
</code><code>  cout &lt;&lt; &quot;Set LHCb Style - Feb 2012&quot; &lt;&lt; endl;
</code><code>  cout &lt;&lt; &quot;-------------------------&quot; &lt;&lt; endl;  
</code><code>  
</code><code>}
</code><code>
</code><code>
</code><code></code></pre>  </div>
</div>
<div id="footer">
  <p style="padding:0; margin:0"><small>Powered by <a href="http://www.websvn.info/">WebSVN</a> 2.3.3 and <a href="http://subversion.tigris.org">Subversion</a> 1.7.14 &nbsp; &nbsp; &#x2713; <a href="http://validator.w3.org/check?uri=https://svnweb.cern.ch/world/wsvn.php?template=%26language=en">XHTML</a> &amp; <a href="http://jigsaw.w3.org/css-validator/validator?uri=https://svnweb.cern.ch/world/wsvn.php?template=%26language=en">CSS</a></small></p>
</div>
</body>
</html>
