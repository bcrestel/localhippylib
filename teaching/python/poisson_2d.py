<!DOCTYPE html>
<html lang="en">
<head>
  
  <meta http-equiv="X-UA-Compatible" content="IE=edge" />
  <meta charset="utf-8">
  <script type="text/javascript">window.NREUM||(NREUM={}),__nr_require=function(e,n,t){function r(t){if(!n[t]){var o=n[t]={exports:{}};e[t][0].call(o.exports,function(n){var o=e[t][1][n];return r(o||n)},o,o.exports)}return n[t].exports}if("function"==typeof __nr_require)return __nr_require;for(var o=0;o<t.length;o++)r(t[o]);return r}({QJf3ax:[function(e,n){function t(e){function n(n,t,a){e&&e(n,t,a),a||(a={});for(var u=c(n),f=u.length,s=i(a,o,r),p=0;f>p;p++)u[p].apply(s,t);return s}function a(e,n){f[e]=c(e).concat(n)}function c(e){return f[e]||[]}function u(){return t(n)}var f={};return{on:a,emit:n,create:u,listeners:c,_events:f}}function r(){return{}}var o="nr@context",i=e("gos");n.exports=t()},{gos:"7eSDFh"}],ee:[function(e,n){n.exports=e("QJf3ax")},{}],3:[function(e,n){function t(e){return function(){r(e,[(new Date).getTime()].concat(i(arguments)))}}var r=e("handle"),o=e(1),i=e(2);"undefined"==typeof window.newrelic&&(newrelic=window.NREUM);var a=["setPageViewName","addPageAction","setCustomAttribute","finished","addToTrace","inlineHit","noticeError"];o(a,function(e,n){window.NREUM[n]=t("api-"+n)}),n.exports=window.NREUM},{1:12,2:13,handle:"D5DuLP"}],gos:[function(e,n){n.exports=e("7eSDFh")},{}],"7eSDFh":[function(e,n){function t(e,n,t){if(r.call(e,n))return e[n];var o=t();if(Object.defineProperty&&Object.keys)try{return Object.defineProperty(e,n,{value:o,writable:!0,enumerable:!1}),o}catch(i){}return e[n]=o,o}var r=Object.prototype.hasOwnProperty;n.exports=t},{}],D5DuLP:[function(e,n){function t(e,n,t){return r.listeners(e).length?r.emit(e,n,t):void(r.q&&(r.q[e]||(r.q[e]=[]),r.q[e].push(n)))}var r=e("ee").create();n.exports=t,t.ee=r,r.q={}},{ee:"QJf3ax"}],handle:[function(e,n){n.exports=e("D5DuLP")},{}],XL7HBI:[function(e,n){function t(e){var n=typeof e;return!e||"object"!==n&&"function"!==n?-1:e===window?0:i(e,o,function(){return r++})}var r=1,o="nr@id",i=e("gos");n.exports=t},{gos:"7eSDFh"}],id:[function(e,n){n.exports=e("XL7HBI")},{}],G9z0Bl:[function(e,n){function t(){var e=d.info=NREUM.info,n=f.getElementsByTagName("script")[0];if(e&&e.licenseKey&&e.applicationID&&n){c(p,function(n,t){n in e||(e[n]=t)});var t="https"===s.split(":")[0]||e.sslForHttp;d.proto=t?"https://":"http://",a("mark",["onload",i()]);var r=f.createElement("script");r.src=d.proto+e.agent,n.parentNode.insertBefore(r,n)}}function r(){"complete"===f.readyState&&o()}function o(){a("mark",["domContent",i()])}function i(){return(new Date).getTime()}var a=e("handle"),c=e(1),u=window,f=u.document;e(2);var s=(""+location).split("?")[0],p={beacon:"bam.nr-data.net",errorBeacon:"bam.nr-data.net",agent:"js-agent.newrelic.com/nr-768.min.js"},d=n.exports={offset:i(),origin:s,features:{}};f.addEventListener?(f.addEventListener("DOMContentLoaded",o,!1),u.addEventListener("load",t,!1)):(f.attachEvent("onreadystatechange",r),u.attachEvent("onload",t)),a("mark",["firstbyte",i()])},{1:12,2:3,handle:"D5DuLP"}],loader:[function(e,n){n.exports=e("G9z0Bl")},{}],12:[function(e,n){function t(e,n){var t=[],o="",i=0;for(o in e)r.call(e,o)&&(t[i]=n(o,e[o]),i+=1);return t}var r=Object.prototype.hasOwnProperty;n.exports=t},{}],13:[function(e,n){function t(e,n,t){n||(n=0),"undefined"==typeof t&&(t=e?e.length:0);for(var r=-1,o=t-n||0,i=Array(0>o?0:o);++r<o;)i[r]=e[n+r];return i}n.exports=t},{}]},{},["G9z0Bl"]);</script>
  <title>
  npetra / math292f15 
  / source  / codes / fenics / forward / poisson_2d.py
 &mdash; Bitbucket
</title>
  


<meta id="bb-canon-url" name="bb-canon-url" content="https://bitbucket.org">

<meta name="description" content=""/>
<meta name="bb-view-name" content="bitbucket.apps.repo2.views.filebrowse">
<meta name="ignore-whitespace" content="False">
<meta name="tab-size" content="None">

<meta name="application-name" content="Bitbucket">
<meta name="apple-mobile-web-app-title" content="Bitbucket">
<meta name="theme-color" content="#205081">
<meta name="msapplication-TileColor" content="#205081">
<meta name="msapplication-TileImage" content="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/img/logos/bitbucket/white-256.png">
<link rel="apple-touch-icon" sizes="192x192" type="image/png" href="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a//img/bitbucket_avatar/192/bitbucket.png">
<link rel="icon" sizes="192x192" type="image/png" href="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a//img/bitbucket_avatar/192/bitbucket.png">
<link rel="icon" sizes="16x16 32x32" type="image/x-icon" href="/favicon.ico">
<link rel="search" type="application/opensearchdescription+xml" href="/opensearch.xml" title="Bitbucket">

  
    
  



  <link rel="stylesheet" href="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/css/entry/vendor.css" />
<link rel="stylesheet" href="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/css/entry/app.css" />


  
  
  
    <link href="/npetra/math292f15/rss?token=9dd30a3d810fba8eb43501abadae84e3" rel="alternate nofollow" type="application/rss+xml" title="RSS feed for math292f15" />
  

</head>
<body class="production aui-page-sidebar aui-sidebar-collapsed"
    data-static-url="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/"
data-base-url="https://bitbucket.org"
data-no-avatar-image="img/default_avatar/16/user_blue.png"
data-current-user="{&quot;username&quot;: &quot;uvilla&quot;, &quot;displayName&quot;: &quot;Umberto Villa&quot;, &quot;uuid&quot;: &quot;{131c19c4-5d93-45d9-9469-091c6bc0f210}&quot;, &quot;firstName&quot;: &quot;Umberto&quot;, &quot;avatarUrl&quot;: &quot;https://bitbucket.org/account/uvilla/avatar/32/?ts=1444892627&quot;, &quot;lastName&quot;: &quot;Villa&quot;, &quot;isTeam&quot;: false, &quot;isSshEnabled&quot;: true, &quot;isKbdShortcutsEnabled&quot;: true, &quot;id&quot;: 5236897, &quot;isAuthenticated&quot;: true}"
data-atlassian-id="{&quot;loginStatusUrl&quot;: &quot;https://id.atlassian.com/profile/rest/profile&quot;}"
data-settings="{&quot;MENTIONS_MIN_QUERY_LENGTH&quot;: 3}"


data-current-repo="{&quot;scm&quot;: &quot;git&quot;, &quot;readOnly&quot;: false, &quot;mainbranch&quot;: {&quot;name&quot;: &quot;master&quot;}, &quot;language&quot;: &quot;&quot;, &quot;owner&quot;: {&quot;username&quot;: &quot;npetra&quot;, &quot;isTeam&quot;: false}, &quot;fullslug&quot;: &quot;npetra/math292f15&quot;, &quot;slug&quot;: &quot;math292f15&quot;, &quot;id&quot;: 13334003, &quot;pygmentsLanguage&quot;: null}"
data-current-cset="909ad9e757141bc9ac07faec7c68a018b9a6a3c9"



>

  
      <script src="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/dist/webpack/early.js"></script>
  


<div id="page">
  <div id="wrapper">
    
  


    <header id="header" role="banner" data-modules="header/tracking">
      
  
  


      <nav class="aui-header aui-dropdown2-trigger-group" role="navigation">
        <div class="aui-header-inner">
          <div class="aui-header-primary">
            
  
  <div class="aui-header-before">
    <button class="app-switcher-trigger aui-dropdown2-trigger aui-dropdown2-trigger-arrowless" aria-controls="app-switcher" aria-haspopup="true" tabindex="0">
      <span class="aui-icon aui-icon-small aui-iconfont-appswitcher">Linked applications</span>
    </button>
    
      <nav id="app-switcher" class="aui-dropdown2 aui-style-default">
        <div class="aui-dropdown2-section blank-slate">
          <h2>Connect Bitbucket with other great Atlassian products:</h2>
          <dl>
            <dt class="jira">JIRA</dt>
            <dd>Project and issue tracking</dd>
            <dt class="confluence">Confluence</dt>
            <dd>Collaboration and content sharing</dd>
            <dt class="bamboo">Bamboo</dt>
            <dd>Continuous integration</dd>
          </dl>
          <ul>
            <li>
              <a href="https://www.atlassian.com/ondemand/signup/?product=jira.ondemand,com.atlassian.bitbucket&utm_source=bitbucket&utm_medium=button&utm_campaign=app_switcher&utm_content=trial_button"
                 class="aui-button aui-button-primary" target="_blank" id="ondemand-trial">Free trial</a>
            </li>
            <li>
              <a href="https://www.atlassian.com/software?utm_source=bitbucket&utm_medium=button&utm_campaign=app_switcher&utm_content=learn_more_button#cloud-products"
                 class="aui-button" target="_blank" id="ondemand-learn-more">Learn more</a>
            </li>
          </ul>
        </div>
      </nav>
    
  </div>


            
              <h1 class="aui-header-logo aui-header-logo-bitbucket " id="logo">
                <a href="/">
                  <span class="aui-header-logo-device">Bitbucket</span>
                </a>
              </h1>
            
            
  
<script id="repo-dropdown-template" type="text/html">
    

[[#hasViewed]]
  <div class="aui-dropdown2-section">
    <strong class="viewed">Recently viewed</strong>
    <ul class="aui-list-truncate">
      [[#viewed]]
        <li class="[[#is_private]]private[[/is_private]][[^is_private]]public[[/is_private]] repository">
          <a href="[[url]]" title="[[owner]]/[[name]]" class="aui-icon-container recently-viewed repo-link">
            <img class="repo-avatar size16" src="[[{avatar}]]" alt="[[owner]]/[[name]] avatar"/>
            [[owner]] / [[name]]
          </a>
        </li>
      [[/viewed]]
    </ul>
  </div>
[[/hasViewed]]
[[#hasUpdated]]
<div class="aui-dropdown2-section">
  <strong class="updated">Recently updated</strong>
  <ul class="aui-list-truncate">
    [[#updated]]
    <li class="[[#is_private]]private[[/is_private]][[^is_private]]public[[/is_private]] repository">
      <a href="[[url]]" title="[[owner]]/[[name]]" class="aui-icon-container recently-updated repo-link">
        <img class="repo-avatar size16" src="[[{avatar}]]" alt="[[owner]]/[[name]] avatar"/>
        [[owner]] / [[name]]
      </a>
    </li>
    [[/updated]]
  </ul>
</div>
[[/hasUpdated]]

  </script>
<script id="snippet-dropdown-template" type="text/html">
    <div class="aui-dropdown2-section">
  <strong>[[sectionTitle]]</strong>
  <ul class="aui-list-truncate">
    [[#snippets]]
      <li>
        <a href="[[links.html.href]]">[[owner.display_name]] / [[name]]</a>
      </li>
    [[/snippets]]
  </ul>
</div>

  </script>
<ul class="aui-nav">
  
    <li>
      <a class="aui-dropdown2-trigger"
         aria-owns="dashboard-dropdown" aria-haspopup="true" href="/dashboard/overview" id="dashboard-dropdown-trigger">
        Dashboard
        <span class="aui-icon-dropdown"></span>
      </a>
      <nav id="dashboard-dropdown" class="aui-dropdown2 aui-style-default">
        <div class="aui-dropdown2-section">
          <ul>
            
              <li>
                <a href="/dashboard/overview">Overview</a>
              </li>
            
              <li>
                <a href="/dashboard/pullrequests">Pull requests</a>
              </li>
            
              <li>
                <a href="/dashboard/issues">Issues</a>
              </li>
            
              <li>
                <a href="/snippets/">Snippets</a>
              </li>
            
          </ul>
        </div>
      </nav>
    </li>
    
      <li>
        <script id="team-dropdown-template" type="text/html">
    

<div class="aui-dropdown2-section primary">
  <ul class="aui-list-truncate">
    [[#teams]]
      <li>
        <a href="/[[username]]/" class="aui-icon-container team-link">
          <span class="aui-avatar aui-avatar-xsmall">
            <span class="aui-avatar-inner">
              <img src="[[avatar]]">
            </span>
          </span>
          [[display_name]]
        </a>
      </li>
    [[/teams]]
  </ul>
</div>

<div class="aui-dropdown2-section">
  <ul>
    <li>
      <a href="/account/create-team/?team_source=header"
          data-modules="registration/create-team-link"
          id="create-team-link">Create team</a>
    </li>
  </ul>
</div>

  </script>
        <a class="aui-dropdown2-trigger" href="#teams-dropdown" id="teams-dropdown-trigger"
          data-modules="header/teams-dropdown"
          aria-owns="teams-dropdown" aria-haspopup="true">
          Teams
          <span class="aui-icon-dropdown"></span>
        </a>
        <div id="teams-dropdown" class="aui-dropdown2 aui-style-default">
          <div class="aui-dropdown2-section blank-slate">
            <p>Organize your team's work and supercharge repository administration.</p>
            <ul>
              <li>
                <a class="aui-button aui-button-primary" id="create-team-blank-slate"
                   href="/account/create-team/?team_source=menu_blank"
                   >Create team</a>
              </li>
            </ul>
          </div>
        </div>
      </li>
    
    <li>
      <a class="aui-dropdown2-trigger" id="repositories-dropdown-trigger"
         aria-owns="repo-dropdown" aria-haspopup="true" href="/repo/mine">
        Repositories
        <span class="aui-icon-dropdown"></span>
      </a>
      <nav id="repo-dropdown" class="aui-dropdown2 aui-style-default">
        <div id="repo-dropdown-recent" data-modules="header/recent-repos"></div>
        <div class="aui-dropdown2-section">
          <ul>
            <li>
              <a href="/repo/create" class="new-repository" id="create-repo-link">
                Create repository
              </a>
            </li>
            <li>
              <a href="/repo/import" class="import-repository" id="import-repo-link">
                Import repository
              </a>
            </li>
          </ul>
        </div>
      </nav>
    </li>
    <li>
      <a class="aui-dropdown2-trigger" id="snippets-dropdown-trigger"
        aria-owns="nav-recent-snippets" aria-haspopup="true" href="/snippets/">
        Snippets
        <span class="aui-icon-dropdown"></span>
      </a>
      <nav id="nav-recent-snippets" class="aui-dropdown2 aui-style-default">
        <div id="snippet-dropdown-recent" class="aui-dropdown2-section"
            data-modules="snippets/recent-list"></div>
        <div class="aui-dropdown2-section">
          <ul>
            <li>
              <a href="/snippets/new">
                Create snippet
              </a>
            </li>
          </ul>
        </div>
      </nav>
    </li>
    <li>
      <button class="aui-button aui-button-primary aui-dropdown2-trigger" aria-owns="create-cta-dropdown" aria-haspopup="true" aria-controls="create-cta-dropdown">
        Create<span class="aui-icon-dropdown"></span>
      </button>
      <nav id="create-cta-dropdown" class="aui-dropdown2 aui-style-default aui-dropdown2-in-header" aria-hidden="true">
        <div class="aui-dropdown2-section">
          <ul>
            
              <li>
                <a href="/repo/create?owner=npetra">Create repository</a>
              </li>
            
              <li>
                <a href="/snippets/new?owner=npetra">Create snippet</a>
              </li>
            
          </ul>
        </div>
      </nav>
    </li>
  
</ul>

          </div>
          <div class="aui-header-secondary">
            
  

<ul role="menu" class="aui-nav">
  
  <li>
    <form action="/repo/all" method="get" class="aui-quicksearch">
      <label for="search-query" class="assistive">owner/repository</label>
      <input id="search-query" class="bb-repo-typeahead" type="text"
             placeholder="Find a repository&hellip;" name="name" autocomplete="off"
             data-bb-typeahead-focus="false">
    </form>
  </li>
  <li id="ace-stp-menu">
    <a id="ace-stp-menu-link" class="aui-nav-link" href="#"
    aria-controls="super-touch-point-dialog"
    data-aui-trigger>
  <span id="ace-stp-menu-icon"
      class="aui-icon aui-icon-small aui-iconfont-help"></span>
</a>
  </li>
  
    
  
  
    <li>
      <a class="aui-dropdown2-trigger aui-dropdown2-trigger-arrowless"
         aria-owns="user-dropdown" aria-haspopup="true" data-container="#header .aui-header-inner"
         href="/uvilla/" title="uvilla" id="user-dropdown-trigger">
        <div class="aui-avatar aui-avatar-small">
            <div class="aui-avatar-inner">
                <img src="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/img/default_avatar/32/user_blue.png" class="deferred-image" data-src-url="https://bitbucket.org/account/uvilla/avatar/32/?ts=1444892627" data-src-url-2x="https://bitbucket.org/account/uvilla/avatar/64/?ts=1444892627" alt="Logged in as uvilla" height="24" width="24">
            </div>
        </div>
      </a>
      <nav id="user-dropdown" class="aui-dropdown2 aui-style-default" aria-hidden="true">
        <div class="aui-dropdown2-section">
          <ul>
            <li>
              <a href="/uvilla/" id="profile-link">View profile</a>
            </li>
            <li>
              <a href="/account/user/uvilla/" id="account-link" data-ct="navbar.dropdown.manage_account">Manage account</a>
            </li>
            <li>
              <a href="/account/user/uvilla/addon-directory" id="account-addons" data-ct="navbar.dropdown.addons">Add-ons</a>
            </li>
            <li>
              <a href="/account/notifications/" id="inbox-link">Inbox</a>
            </li>
            <li>
              <a href="#general-invite" class="general-invite-link" id="general-invite-dropdown" data-modules="header/general-invite">Invite a friend</a>
              <script id="general-invite-template" type="text/html">
    
<div id="general-invite">
  <form class="aui invitation-form" id="general-invite-form" method="post" novalidate>
    <div class="error aui-message hidden">
      <span class="aui-icon icon-error"></span>
      <div class="message"></div>
    </div>
    <div class="field-group">
      <label for="id_general_email_address">Email address</label>
      <input type="email" id="id_general_email_address" name="email-address" class="text long-field">
    </div>
  </form>
</div>

  </script>
            </li>
          </ul>
        </div>
        <div class="aui-dropdown2-section">
          <ul>
            <li>
              
                <a href="/account/signout/" id="log-out-link">Log out</a>
              
            </li>
          </ul>
        </div>
      </nav>
    </li>
  
</ul>

          </div>
        </div>
      </nav>
    </header>

    
  

<header id="account-warning" role="banner" data-modules="header/account-warning"
        class="aui-message-banner warning
        ">
  <div class="aui-message-banner-inner">
    <span class="aui-icon aui-icon-warning"></span>
    <span class="message">
    
    </span>
  </div>
</header>



    
  
<header id="aui-message-bar">
  
</header>


    <div id="content" role="main">
      
  

<div class="aui-sidebar repo-sidebar"
    data-modules="components/sidebar,components/repo-sidebar"
   aria-expanded="false">
  <div class="aui-sidebar-wrapper">
    <div class="aui-sidebar-body">
      <header class="aui-page-header">
        <div class="aui-page-header-inner">
          <div class="aui-page-header-image">
            <a href="/npetra/math292f15" id="repo-avatar-link" class="repo-link">
              <span class="aui-avatar aui-avatar-large aui-avatar-project">
                <span class="aui-avatar-inner">
                  <img  src="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/img/language-avatars/default_64.png" class="deferred-image" data-src-url="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/img/language-avatars/default_64.png" data-src-url-2x="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/img/language-avatars/default_128.png" alt="">
                </span>
              </span>
            </a>
          </div>
          <div class="aui-page-header-main">
            <ol class="aui-nav aui-nav-breadcrumbs">
              <li>
                <a href="/npetra/" id="repo-owner-link">npetra</a>
              </li>
            </ol>
            <h1>
              
                <span class="aui-icon aui-icon-small aui-iconfont-locked"></span>
              
              <a href="/npetra/math292f15" title="math292f15" class="entity-name">math292f15</a>
            </h1>
          </div>
        </div>
      </header>
      <nav class="aui-navgroup aui-navgroup-vertical">
        <div class="aui-navgroup-inner">
          
            
              <div class="aui-sidebar-group aui-sidebar-group-actions repository-actions forks-enabled can-create">
                <div class="aui-nav-heading">
                  <strong>Actions</strong>
                </div>
                <ul id="repo-actions" class="aui-nav">
                  
                  
                    <li>
                      <a id="repo-clone-button" class="aui-nav-item "
                        href="#clone"
                        data-ct="sidebar.actions.repository.clone"
                        data-ct-data=""
                        data-modules="components/clone/clone-dialog"
                        target="_self">
                        
                          <span class="aui-icon aui-icon-large icon-clone"></span>
                        
                        <span class="aui-nav-item-label">Clone</span>
                      </a>
                    </li>
                  
                    <li>
                      <a id="repo-create-branch-link" class="aui-nav-item create-branch-button"
                        href="/npetra/math292f15/branch"
                        data-ct="sidebar.actions.repository.create_branch"
                        data-ct-data=""
                        data-modules="repo/create-branch"
                        target="_self">
                        
                          <span class="aui-icon aui-icon-large icon-create-branch"></span>
                        
                        <span class="aui-nav-item-label">Create branch</span>
                      </a>
                    </li>
                  
                    <li>
                      <a id="repo-create-pull-request-link" class="aui-nav-item "
                        href="/npetra/math292f15/pull-requests/new"
                        data-ct="sidebar.actions.create_pullrequest"
                        data-ct-data=""
                        
                        target="_self">
                        
                          <span class="aui-icon aui-icon-large icon-create-pull-request"></span>
                        
                        <span class="aui-nav-item-label">Create pull request</span>
                      </a>
                    </li>
                  
                    <li>
                      <a id="repo-compare-link" class="aui-nav-item "
                        href="/npetra/math292f15/branches/compare"
                        data-ct="sidebar.actions.repository.compare"
                        data-ct-data=""
                        
                        target="_self">
                        
                          <span class="aui-icon aui-icon-large aui-icon-small aui-iconfont-devtools-compare"></span>
                        
                        <span class="aui-nav-item-label">Compare</span>
                      </a>
                    </li>
                  
                    <li>
                      <a id="repo-fork-link" class="aui-nav-item "
                        href="/npetra/math292f15/fork"
                        data-ct="sidebar.actions.repository.fork"
                        data-ct-data=""
                        
                        target="_self">
                        
                          <span class="aui-icon aui-icon-large icon-fork"></span>
                        
                        <span class="aui-nav-item-label">Fork</span>
                      </a>
                    </li>
                  
                </ul>
              </div>
          
          <div class="aui-sidebar-group aui-sidebar-group-tier-one repository-sections">
            <div class="aui-nav-heading">
              <strong>Navigation</strong>
            </div>
            <ul class="aui-nav">
              
              
                <li>
                  <a id="repo-overview-link" class="aui-nav-item "
                    href="/npetra/math292f15/overview"
                    data-ct="sidebar.navigation.repository.overview"
                    data-ct-data=""
                    
                    target="_self"
                    >
                    
                    <span class="aui-icon aui-icon-large icon-overview"></span>
                    <span class="aui-nav-item-label">Overview</span>
                  </a>
                </li>
              
                <li class="aui-nav-selected">
                  <a id="repo-source-link" class="aui-nav-item "
                    href="/npetra/math292f15/src"
                    data-ct="sidebar.navigation.repository.source"
                    data-ct-data=""
                    
                    target="_self"
                    >
                    
                    <span class="aui-icon aui-icon-large icon-source"></span>
                    <span class="aui-nav-item-label">Source</span>
                  </a>
                </li>
              
                <li>
                  <a id="repo-commits-link" class="aui-nav-item "
                    href="/npetra/math292f15/commits/"
                    data-ct="sidebar.navigation.repository.commits"
                    data-ct-data=""
                    
                    target="_self"
                    >
                    
                    <span class="aui-icon aui-icon-large icon-commits"></span>
                    <span class="aui-nav-item-label">Commits</span>
                  </a>
                </li>
              
                <li>
                  <a id="repo-branches-link" class="aui-nav-item "
                    href="/npetra/math292f15/branches/"
                    data-ct="sidebar.navigation.repository.branches"
                    data-ct-data=""
                    
                    target="_self"
                    >
                    
                    <span class="aui-icon aui-icon-large icon-branches"></span>
                    <span class="aui-nav-item-label">Branches</span>
                  </a>
                </li>
              
                <li>
                  <a id="repo-pullrequests-link" class="aui-nav-item "
                    href="/npetra/math292f15/pull-requests/"
                    data-ct="sidebar.navigation.repository.pullrequests"
                    data-ct-data=""
                    
                    target="_self"
                    >
                    
                    <span class="aui-icon aui-icon-large icon-pull-requests"></span>
                    <span class="aui-nav-item-label">Pull requests</span>
                  </a>
                </li>
              
                <li>
                  <a id="repo-downloads-link" class="aui-nav-item "
                    href="/npetra/math292f15/downloads"
                    data-ct="sidebar.navigation.repository.downloads"
                    data-ct-data=""
                    
                    target="_self"
                    >
                    
                    <span class="aui-icon aui-icon-large icon-downloads"></span>
                    <span class="aui-nav-item-label">Downloads</span>
                  </a>
                </li>
              
            </ul>
          </div>
          <div class="aui-sidebar-group aui-sidebar-group-tier-one repository-settings">
            <div class="aui-nav-heading">
              <strong class="assistive">Settings</strong>
            </div>
            <ul class="aui-nav">
              
              
            </ul>
          </div>
          
        </div>
      </nav>
    </div>
    <div class="aui-sidebar-footer">
      <a class="aui-sidebar-toggle aui-sidebar-footer-tipsy aui-button aui-button-subtle"><span class="aui-icon"></span></a>
    </div>
  </div>
  
    <script id="share-repo-template" type="text/html">
    

<div class="clearfix">
   <span class="repo-avatar-container size32" title="npetra/math292f15">
  <img alt="npetra/math292f15" src="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/img/language-avatars/default_32.png" class="deferred-image" data-src-url="https://bitbucket.org/npetra/math292f15/avatar/32/?ts=1446154377" data-src-url-2x="https://bitbucket.org/npetra/math292f15/avatar/64/?ts=1446154377">
</span>

  <span class="repo-name-container">
    npetra / math292f15
  </span>
</div>
<p class="helptext">
  
    Existing users are granted access to this repository immediately.
    New users will be sent an invitation.
  
</p>
<div class="manage-repo-link"
  data-manage-link="/npetra/math292f15/admin/access"></div>
<div class="share-form"></div>

  </script>
    <script id="share-dialog-template" type="text/html">
    <div class="aui-tabs horizontal-tabs">
  <ul class="tabs-menu">
    [[#panels]]
      <li class="menu-item">
        <a href="#[[tabId]]"><strong>[[display]]</strong></a>
      </li>
    [[/panels]]
  </ul>
  [[#panels]]
    <div class="tabs-pane" id="[[tabId]]"></div>
  [[/panels]]
</div>

  </script>
  
  <div id="repo-clone-dialog" class="clone-dialog hidden">
  
  
  <div class="clone-url" data-modules="components/clone/url-dropdown">
  <div class="aui-buttons">
    <a href="https://uvilla@bitbucket.org/npetra/math292f15.git"
      class="aui-button aui-dropdown2-trigger" aria-haspopup="true"
      aria-owns="clone-url-dropdown-header">
      <span class="dropdown-text">HTTPS</span>
    </a>
    <div id="clone-url-dropdown-header"
        class="clone-url-dropdown aui-dropdown2 aui-style-default"
        data-aui-alignment="bottom left">
      <ul class="aui-list-truncate">
        <li>
          <a href="https://uvilla@bitbucket.org/npetra/math292f15.git"
            
              data-command="git clone https://uvilla@bitbucket.org/npetra/math292f15.git"
            
            class="item-link https">HTTPS
          </a>
        </li>
        <li>
          <a href="ssh://git@bitbucket.org/npetra/math292f15.git"
            
              data-command="git clone git@bitbucket.org:npetra/math292f15.git"
            
            class="item-link ssh">SSH
          </a>
        </li>
      </ul>
    </div>
    <input type="text" readonly="readonly" class="clone-url-input"
      value="git clone https://uvilla@bitbucket.org/npetra/math292f15.git">
  </div>
  
  <p>Need help cloning? Learn how to
       <a href="https://confluence.atlassian.com/x/4whODQ" target="_blank">clone a repository</a>.</p>
  
</div>
  
  <div class="sourcetree-callout clone-in-sourcetree"
  data-modules="components/clone/clone-in-sourcetree"
  data-https-url="https://uvilla@bitbucket.org/npetra/math292f15.git"
  data-ssh-url="ssh://git@bitbucket.org/npetra/math292f15.git">

  <div>
    <button class="aui-button aui-button-primary">
      
        Clone in SourceTree
      
    </button>
  </div>

  <p class="windows-text">
    
      <a href="http://www.sourcetreeapp.com/?utm_source=internal&amp;utm_medium=link&amp;utm_campaign=clone_repo_win" target="_blank">Atlassian SourceTree</a>
      is a free Git and Mercurial client for Windows.
    
  </p>
  <p class="mac-text">
    
      <a href="http://www.sourcetreeapp.com/?utm_source=internal&amp;utm_medium=link&amp;utm_campaign=clone_repo_mac" target="_blank">Atlassian SourceTree</a>
      is a free Git and Mercurial client for Mac.
    
  </p>
</div>
</div>
</div>

      
  <div class="aui-page-panel ">
    



    <div class="aui-page-panel-inner">
      <div id="repo-content" class="aui-page-panel-content" data-modules="repo/index">
        <div class="aui-group repo-page-header">
          <div class="aui-item section-title">
            <h1>Source</h1>
          </div>
          <div class="aui-item page-actions">
            
          </div>
        </div>
        
  <div id="source-container" class="maskable" data-modules="repo/source/index">
    



<header id="source-path">
  <div class="labels labels-csv">
    
      <div class="aui-buttons">
        <button data-branches-tags-url="/api/1.0/repositories/npetra/math292f15/branches-tags"
                data-modules="components/branch-dialog"
                class="aui-button branch-dialog-trigger" title="master">
          
            
              <span class="aui-icon aui-icon-small aui-iconfont-devtools-branch">Branch</span>
            
            <span class="name">master</span>
          
          <span class="aui-icon-dropdown"></span>
        </button>
        <button class="aui-button" id="checkout-branch-button"
                title="Check out this branch">
          <span class="aui-icon aui-icon-small aui-iconfont-devtools-clone">Check out branch</span>
          <span class="aui-icon-dropdown"></span>
        </button>
      </div>
      <script id="branch-checkout-template" type="text/html">
  

<div id="checkout-branch-contents">
  <div class="command-line">
    <p>
      Check out this branch on your local machine to begin working on it.
    </p>
    <input type="text" class="checkout-command" readonly="readonly"
        
           value="git fetch && git checkout [[branchName]]"
        
        >
  </div>
  
    <div class="sourcetree-callout clone-in-sourcetree"
  data-modules="components/clone/clone-in-sourcetree"
  data-https-url="https://uvilla@bitbucket.org/npetra/math292f15.git"
  data-ssh-url="ssh://git@bitbucket.org/npetra/math292f15.git">

  <div>
    <button class="aui-button aui-button-primary">
      
        Check out in SourceTree
      
    </button>
  </div>

  <p class="windows-text">
    
      <a href="http://www.sourcetreeapp.com/?utm_source=internal&amp;utm_medium=link&amp;utm_campaign=clone_repo_win" target="_blank">Atlassian SourceTree</a>
      is a free Git and Mercurial client for Windows.
    
  </p>
  <p class="mac-text">
    
      <a href="http://www.sourcetreeapp.com/?utm_source=internal&amp;utm_medium=link&amp;utm_campaign=clone_repo_mac" target="_blank">Atlassian SourceTree</a>
      is a free Git and Mercurial client for Mac.
    
  </p>
</div>
  
</div>

</script>
    
  </div>
  <div class="secondary-actions">
    <div class="aui-buttons">
      
        <a href="/npetra/math292f15/src/909ad9e75714/codes/fenics/forward/poisson_2d.py?at=master"
           class="aui-button pjax-trigger" aria-pressed="true">
          Source
        </a>
        <a href="/npetra/math292f15/diff/codes/fenics/forward/poisson_2d.py?diff2=909ad9e75714&at=master"
           class="aui-button pjax-trigger"
           title="Diff to previous change">
          Diff
        </a>
        <a href="/npetra/math292f15/history-node/909ad9e75714/codes/fenics/forward/poisson_2d.py?at=master"
           class="aui-button pjax-trigger">
          History
        </a>
      
    </div>
  </div>
  <h1>
    
      
        <a href="/npetra/math292f15/src/909ad9e75714?at=master"
          class="pjax-trigger root" title="npetra/math292f15 at 909ad9e75714">math292f15</a> /
      
      
        
          
            
              <a href="/npetra/math292f15/src/909ad9e75714/codes/?at=master"
                class="pjax-trigger directory-name">codes</a> /
            
          
        
      
        
          
            
              <a href="/npetra/math292f15/src/909ad9e75714/codes/fenics/?at=master"
                class="pjax-trigger directory-name">fenics</a> /
            
          
        
      
        
          
            
              <a href="/npetra/math292f15/src/909ad9e75714/codes/fenics/forward/?at=master"
                class="pjax-trigger directory-name">forward</a> /
            
          
        
      
        
          
            
              <span class="file-name">poisson_2d.py</span>
            
          
        
      
    
  </h1>
  
    
    
  
  <div class="clearfix"></div>
</header>


    
      
    

    <div id="editor-container" class="maskable"
         data-modules="repo/source/editor"
         data-owner="npetra"
         data-slug="math292f15"
         data-is-writer="true"
         data-has-push-access="true"
         data-hash="909ad9e757141bc9ac07faec7c68a018b9a6a3c9"
         data-branch="master"
         data-path="codes/fenics/forward/poisson_2d.py"
         data-source-url="/api/internal//npetra/math292f15/src/909ad9e757141bc9ac07faec7c68a018b9a6a3c9/codes/fenics/forward/poisson_2d.py">
      <div id="source-view" class="file-source-container" data-modules="repo/source/view-file">
        <div class="toolbar">
          <div class="primary">
            <div class="aui-buttons">
              
                <button id="file-history-trigger" class="aui-button aui-button-light changeset-info"
                        data-changeset="909ad9e757141bc9ac07faec7c68a018b9a6a3c9"
                        data-path="codes/fenics/forward/poisson_2d.py"
                        data-current="909ad9e757141bc9ac07faec7c68a018b9a6a3c9">
                  
                     

  <div class="aui-avatar aui-avatar-xsmall">
    <div class="aui-avatar-inner">
      <img src="https://bitbucket.org/account/npetra/avatar/16/?ts=0">
    </div>
  </div>
  <span class="changeset-hash">909ad9e</span>
  <time datetime="2015-10-29T21:32:56+00:00" class="timestamp"></time>
  <span class="aui-icon-dropdown"></span>

                  
                </button>
              
            </div>
          <a href="/npetra/math292f15/full-commit/909ad9e75714/codes/fenics/forward/poisson_2d.py" id="full-commit-link"
              title="View full commit 909ad9e">Full commit</a>
          </div>
            <div class="secondary">
              
              <div class="aui-buttons">
                
                  <a href="/npetra/math292f15/annotate/909ad9e757141bc9ac07faec7c68a018b9a6a3c9/codes/fenics/forward/poisson_2d.py?at=master"
                      class="aui-button aui-button-light pjax-trigger">Blame</a>
                
                
                <a href="/npetra/math292f15/raw/909ad9e757141bc9ac07faec7c68a018b9a6a3c9/codes/fenics/forward/poisson_2d.py" class="aui-button aui-button-light">Raw</a>
              </div>
              
                <div class="aui-buttons">
                  <button class="edit-button aui-button aui-button-light aui-button-split-main">
                    Edit
                    
                  </button>
                  <button class="aui-button aui-button-light aui-dropdown2-trigger aui-button-split-more" aria-owns="split-container-dropdown" aria-haspopup="true" >More file actions</button>
                </div>
                <div id="split-container-dropdown" class="aui-dropdown2 aui-style-default" data-container="#editor-container">
                  <ul class="aui-list-truncate">
                    <li><a href="#" data-modules="repo/source/rename-file" class="rename-link">Rename</a></li>
                    <li><a href="#" data-modules="repo/source/delete-file" class="delete-link">Delete</a></li>
                  </ul>
                </div>
              
            </div>

            <div id="fileview-dropdown"
                class="aui-dropdown2 aui-style-default"
                data-fileview-container="#fileview-container"
                
                
                data-fileview-button="#fileview-trigger"
                data-modules="connect/fileview">
              <div class="aui-dropdown2-section">
                <ul>
                  <li>
                    <a class="aui-dropdown2-radio aui-dropdown2-checked"
                        data-fileview-id="-1"
                        data-fileview-loaded="true"
                        data-fileview-target="#fileview-original"
                        data-fileview-connection-key=""
                        data-fileview-module-key="file-view-default">
                      Default File Viewer
                    </a>
                  </li>
                  
                </ul>
              </div>
            </div>

          <div class="clearfix"></div>
        </div>
        <div id="fileview-container">
          <div id="fileview-original"
              class="fileview"
              data-modules="repo/source/highlight-lines"
              data-fileview-loaded="true">
            


  
    <div class="file-source">
      <table class="codehilite highlighttable"><tr><td class="linenos"><div class="linenodiv"><pre><a href="#poisson_2d.py-1">  1</a>
<a href="#poisson_2d.py-2">  2</a>
<a href="#poisson_2d.py-3">  3</a>
<a href="#poisson_2d.py-4">  4</a>
<a href="#poisson_2d.py-5">  5</a>
<a href="#poisson_2d.py-6">  6</a>
<a href="#poisson_2d.py-7">  7</a>
<a href="#poisson_2d.py-8">  8</a>
<a href="#poisson_2d.py-9">  9</a>
<a href="#poisson_2d.py-10"> 10</a>
<a href="#poisson_2d.py-11"> 11</a>
<a href="#poisson_2d.py-12"> 12</a>
<a href="#poisson_2d.py-13"> 13</a>
<a href="#poisson_2d.py-14"> 14</a>
<a href="#poisson_2d.py-15"> 15</a>
<a href="#poisson_2d.py-16"> 16</a>
<a href="#poisson_2d.py-17"> 17</a>
<a href="#poisson_2d.py-18"> 18</a>
<a href="#poisson_2d.py-19"> 19</a>
<a href="#poisson_2d.py-20"> 20</a>
<a href="#poisson_2d.py-21"> 21</a>
<a href="#poisson_2d.py-22"> 22</a>
<a href="#poisson_2d.py-23"> 23</a>
<a href="#poisson_2d.py-24"> 24</a>
<a href="#poisson_2d.py-25"> 25</a>
<a href="#poisson_2d.py-26"> 26</a>
<a href="#poisson_2d.py-27"> 27</a>
<a href="#poisson_2d.py-28"> 28</a>
<a href="#poisson_2d.py-29"> 29</a>
<a href="#poisson_2d.py-30"> 30</a>
<a href="#poisson_2d.py-31"> 31</a>
<a href="#poisson_2d.py-32"> 32</a>
<a href="#poisson_2d.py-33"> 33</a>
<a href="#poisson_2d.py-34"> 34</a>
<a href="#poisson_2d.py-35"> 35</a>
<a href="#poisson_2d.py-36"> 36</a>
<a href="#poisson_2d.py-37"> 37</a>
<a href="#poisson_2d.py-38"> 38</a>
<a href="#poisson_2d.py-39"> 39</a>
<a href="#poisson_2d.py-40"> 40</a>
<a href="#poisson_2d.py-41"> 41</a>
<a href="#poisson_2d.py-42"> 42</a>
<a href="#poisson_2d.py-43"> 43</a>
<a href="#poisson_2d.py-44"> 44</a>
<a href="#poisson_2d.py-45"> 45</a>
<a href="#poisson_2d.py-46"> 46</a>
<a href="#poisson_2d.py-47"> 47</a>
<a href="#poisson_2d.py-48"> 48</a>
<a href="#poisson_2d.py-49"> 49</a>
<a href="#poisson_2d.py-50"> 50</a>
<a href="#poisson_2d.py-51"> 51</a>
<a href="#poisson_2d.py-52"> 52</a>
<a href="#poisson_2d.py-53"> 53</a>
<a href="#poisson_2d.py-54"> 54</a>
<a href="#poisson_2d.py-55"> 55</a>
<a href="#poisson_2d.py-56"> 56</a>
<a href="#poisson_2d.py-57"> 57</a>
<a href="#poisson_2d.py-58"> 58</a>
<a href="#poisson_2d.py-59"> 59</a>
<a href="#poisson_2d.py-60"> 60</a>
<a href="#poisson_2d.py-61"> 61</a>
<a href="#poisson_2d.py-62"> 62</a>
<a href="#poisson_2d.py-63"> 63</a>
<a href="#poisson_2d.py-64"> 64</a>
<a href="#poisson_2d.py-65"> 65</a>
<a href="#poisson_2d.py-66"> 66</a>
<a href="#poisson_2d.py-67"> 67</a>
<a href="#poisson_2d.py-68"> 68</a>
<a href="#poisson_2d.py-69"> 69</a>
<a href="#poisson_2d.py-70"> 70</a>
<a href="#poisson_2d.py-71"> 71</a>
<a href="#poisson_2d.py-72"> 72</a>
<a href="#poisson_2d.py-73"> 73</a>
<a href="#poisson_2d.py-74"> 74</a>
<a href="#poisson_2d.py-75"> 75</a>
<a href="#poisson_2d.py-76"> 76</a>
<a href="#poisson_2d.py-77"> 77</a>
<a href="#poisson_2d.py-78"> 78</a>
<a href="#poisson_2d.py-79"> 79</a>
<a href="#poisson_2d.py-80"> 80</a>
<a href="#poisson_2d.py-81"> 81</a>
<a href="#poisson_2d.py-82"> 82</a>
<a href="#poisson_2d.py-83"> 83</a>
<a href="#poisson_2d.py-84"> 84</a>
<a href="#poisson_2d.py-85"> 85</a>
<a href="#poisson_2d.py-86"> 86</a>
<a href="#poisson_2d.py-87"> 87</a>
<a href="#poisson_2d.py-88"> 88</a>
<a href="#poisson_2d.py-89"> 89</a>
<a href="#poisson_2d.py-90"> 90</a>
<a href="#poisson_2d.py-91"> 91</a>
<a href="#poisson_2d.py-92"> 92</a>
<a href="#poisson_2d.py-93"> 93</a>
<a href="#poisson_2d.py-94"> 94</a>
<a href="#poisson_2d.py-95"> 95</a>
<a href="#poisson_2d.py-96"> 96</a>
<a href="#poisson_2d.py-97"> 97</a>
<a href="#poisson_2d.py-98"> 98</a>
<a href="#poisson_2d.py-99"> 99</a>
<a href="#poisson_2d.py-100">100</a>
<a href="#poisson_2d.py-101">101</a>
<a href="#poisson_2d.py-102">102</a>
<a href="#poisson_2d.py-103">103</a>
<a href="#poisson_2d.py-104">104</a>
<a href="#poisson_2d.py-105">105</a>
<a href="#poisson_2d.py-106">106</a>
<a href="#poisson_2d.py-107">107</a>
<a href="#poisson_2d.py-108">108</a>
<a href="#poisson_2d.py-109">109</a>
<a href="#poisson_2d.py-110">110</a>
<a href="#poisson_2d.py-111">111</a>
<a href="#poisson_2d.py-112">112</a>
<a href="#poisson_2d.py-113">113</a>
<a href="#poisson_2d.py-114">114</a>
<a href="#poisson_2d.py-115">115</a>
<a href="#poisson_2d.py-116">116</a>
<a href="#poisson_2d.py-117">117</a>
<a href="#poisson_2d.py-118">118</a>
<a href="#poisson_2d.py-119">119</a>
<a href="#poisson_2d.py-120">120</a>
<a href="#poisson_2d.py-121">121</a>
<a href="#poisson_2d.py-122">122</a>
<a href="#poisson_2d.py-123">123</a>
<a href="#poisson_2d.py-124">124</a>
<a href="#poisson_2d.py-125">125</a>
<a href="#poisson_2d.py-126">126</a></pre></div></td><td class="code"><div class="codehilite highlight"><pre><a name="poisson_2d.py-1"></a><span class="sd">&quot;&quot;&quot;</span>
<a name="poisson_2d.py-2"></a><span class="sd">Solve the boundary value problem (BVP) that has the form:</span>
<a name="poisson_2d.py-3"></a>
<a name="poisson_2d.py-4"></a><span class="sd">- div (k grad u) = f,      in Omega,</span>
<a name="poisson_2d.py-5"></a><span class="sd">               u = u0,     on Gamma_D = Gamma_left U Gamma_right</span>
<a name="poisson_2d.py-6"></a><span class="sd">            du/dn = sigma, on Gamma_N = Gamma_top U Gamma_bottom</span>
<a name="poisson_2d.py-7"></a>
<a name="poisson_2d.py-8"></a><span class="sd">where Omega = (0,1) x (0,1), Gamma_D and and Omega_N are the union of</span>
<a name="poisson_2d.py-9"></a><span class="sd">the left and right, and top and bottom boundaries of Omega,</span>
<a name="poisson_2d.py-10"></a><span class="sd">respectively. The functions k(x,y) = 1,</span>
<a name="poisson_2d.py-11"></a><span class="sd">f(x,y)     = (4.0*pi*pi+pi*pi/4.0)*(sin(2*pi*x[0])*sin((pi/2.0)*x[1])</span>
<a name="poisson_2d.py-12"></a><span class="sd">u0(y)      = 0, on Gamma_D</span>
<a name="poisson_2d.py-13"></a><span class="sd">sigma(x)   = -pi/2(sin(2*pi*x[0])) on Gamma_bottom,</span>
<a name="poisson_2d.py-14"></a><span class="sd">           = 0,                    on Gamma_top.</span>
<a name="poisson_2d.py-15"></a>
<a name="poisson_2d.py-16"></a><span class="sd">The exact solution is</span>
<a name="poisson_2d.py-17"></a><span class="sd">u_e(x,y) = sin(2*pi*x[0])*sin((pi/2.0)*x[1])</span>
<a name="poisson_2d.py-18"></a>
<a name="poisson_2d.py-19"></a><span class="sd">This example is based on the codes provided in the FEniCS examples, see</span>
<a name="poisson_2d.py-20"></a><span class="sd">http://fenicsproject.org/documentation/tutorial/index.html</span>
<a name="poisson_2d.py-21"></a><span class="sd">Last update: October 13, 2015</span>
<a name="poisson_2d.py-22"></a><span class="sd">Author: Noemi Petra (npetra@ucmerced.edu)</span>
<a name="poisson_2d.py-23"></a><span class="sd">&quot;&quot;&quot;</span>
<a name="poisson_2d.py-24"></a><span class="kn">from</span> <span class="nn">dolfin</span> <span class="kn">import</span> <span class="o">*</span>
<a name="poisson_2d.py-25"></a><span class="kn">import</span> <span class="nn">numpy</span>
<a name="poisson_2d.py-26"></a><span class="kn">import</span> <span class="nn">sys</span>
<a name="poisson_2d.py-27"></a>
<a name="poisson_2d.py-28"></a><span class="n">set_log_active</span><span class="p">(</span><span class="bp">False</span><span class="p">)</span>
<a name="poisson_2d.py-29"></a><span class="k">print</span> <span class="s">&quot;Messages and Warnings Disabled&quot;</span>
<a name="poisson_2d.py-30"></a>
<a name="poisson_2d.py-31"></a><span class="n">nx</span> <span class="o">=</span> <span class="nb">int</span><span class="p">(</span><span class="n">sys</span><span class="o">.</span><span class="n">argv</span><span class="p">[</span><span class="mi">1</span><span class="p">])</span>
<a name="poisson_2d.py-32"></a><span class="n">ny</span> <span class="o">=</span> <span class="n">nx</span>
<a name="poisson_2d.py-33"></a><span class="n">fem_degree</span> <span class="o">=</span> <span class="mi">1</span>
<a name="poisson_2d.py-34"></a>
<a name="poisson_2d.py-35"></a><span class="c"># Create mesh and define function space</span>
<a name="poisson_2d.py-36"></a><span class="n">mesh</span> <span class="o">=</span> <span class="n">RectangleMesh</span><span class="p">(</span><span class="mi">0</span><span class="p">,</span> <span class="mi">0</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="n">nx</span><span class="p">,</span> <span class="n">nx</span><span class="p">)</span>
<a name="poisson_2d.py-37"></a><span class="n">V</span>  <span class="o">=</span> <span class="n">FunctionSpace</span><span class="p">(</span><span class="n">mesh</span><span class="p">,</span> <span class="s">&#39;Lagrange&#39;</span><span class="p">,</span> <span class="n">degree</span><span class="o">=</span><span class="n">fem_degree</span><span class="p">)</span>
<a name="poisson_2d.py-38"></a>
<a name="poisson_2d.py-39"></a><span class="c"># Create mesh function over cell facets</span>
<a name="poisson_2d.py-40"></a><span class="n">boundary_parts</span> <span class="o">=</span> <span class="n">MeshFunction</span><span class="p">(</span><span class="s">&quot;size_t&quot;</span><span class="p">,</span> <span class="n">mesh</span><span class="p">,</span> <span class="n">mesh</span><span class="o">.</span><span class="n">topology</span><span class="p">()</span><span class="o">.</span><span class="n">dim</span><span class="p">()</span><span class="o">-</span><span class="mi">1</span><span class="p">)</span>
<a name="poisson_2d.py-41"></a><span class="n">boundary_parts</span><span class="o">.</span><span class="n">set_all</span><span class="p">(</span><span class="mi">4</span><span class="p">)</span>
<a name="poisson_2d.py-42"></a>
<a name="poisson_2d.py-43"></a><span class="c"># Mark lower boundary facets as subdomain 0</span>
<a name="poisson_2d.py-44"></a><span class="k">class</span> <span class="nc">BottomBoundary</span><span class="p">(</span><span class="n">SubDomain</span><span class="p">):</span>
<a name="poisson_2d.py-45"></a>    <span class="k">def</span> <span class="nf">inside</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span> <span class="n">x</span><span class="p">,</span> <span class="n">on_boundary</span><span class="p">):</span>
<a name="poisson_2d.py-46"></a>        <span class="n">tol</span> <span class="o">=</span> <span class="mf">1E-16</span>   <span class="c"># tolerance for coordinate comparisons</span>
<a name="poisson_2d.py-47"></a>        <span class="k">return</span> <span class="n">on_boundary</span> <span class="ow">and</span> <span class="nb">abs</span><span class="p">(</span><span class="n">x</span><span class="p">[</span><span class="mi">1</span><span class="p">])</span> <span class="o">&lt;</span> <span class="n">tol</span>
<a name="poisson_2d.py-48"></a><span class="n">Gamma_bottom</span> <span class="o">=</span> <span class="n">BottomBoundary</span><span class="p">()</span>
<a name="poisson_2d.py-49"></a><span class="n">Gamma_bottom</span><span class="o">.</span><span class="n">mark</span><span class="p">(</span><span class="n">boundary_parts</span><span class="p">,</span> <span class="mi">0</span><span class="p">)</span>
<a name="poisson_2d.py-50"></a>
<a name="poisson_2d.py-51"></a><span class="c"># Mark upper boundary facets as subdomain 1</span>
<a name="poisson_2d.py-52"></a><span class="k">class</span> <span class="nc">TopBoundary</span><span class="p">(</span><span class="n">SubDomain</span><span class="p">):</span>
<a name="poisson_2d.py-53"></a>    <span class="k">def</span> <span class="nf">inside</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span> <span class="n">x</span><span class="p">,</span> <span class="n">on_boundary</span><span class="p">):</span>
<a name="poisson_2d.py-54"></a>        <span class="n">tol</span> <span class="o">=</span> <span class="mf">1E-16</span>   <span class="c"># tolerance for coordinate comparisons</span>
<a name="poisson_2d.py-55"></a>        <span class="k">return</span> <span class="n">on_boundary</span> <span class="ow">and</span> <span class="nb">abs</span><span class="p">(</span><span class="n">x</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span> <span class="o">-</span> <span class="mi">1</span><span class="p">)</span> <span class="o">&lt;</span> <span class="n">tol</span>
<a name="poisson_2d.py-56"></a><span class="n">Gamma_top</span> <span class="o">=</span> <span class="n">TopBoundary</span><span class="p">()</span>
<a name="poisson_2d.py-57"></a><span class="n">Gamma_top</span><span class="o">.</span><span class="n">mark</span><span class="p">(</span><span class="n">boundary_parts</span><span class="p">,</span> <span class="mi">1</span><span class="p">)</span>
<a name="poisson_2d.py-58"></a>
<a name="poisson_2d.py-59"></a><span class="c"># Mark left boundary as subdomain 2</span>
<a name="poisson_2d.py-60"></a><span class="k">class</span> <span class="nc">LeftBoundary</span><span class="p">(</span><span class="n">SubDomain</span><span class="p">):</span>
<a name="poisson_2d.py-61"></a>    <span class="k">def</span> <span class="nf">inside</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span> <span class="n">x</span><span class="p">,</span> <span class="n">on_boundary</span><span class="p">):</span>
<a name="poisson_2d.py-62"></a>        <span class="n">tol</span> <span class="o">=</span> <span class="mf">1E-16</span>   <span class="c"># tolerance for coordinate comparisons</span>
<a name="poisson_2d.py-63"></a>        <span class="k">return</span> <span class="n">on_boundary</span> <span class="ow">and</span> <span class="nb">abs</span><span class="p">(</span><span class="n">x</span><span class="p">[</span><span class="mi">0</span><span class="p">])</span> <span class="o">&lt;</span> <span class="n">tol</span>
<a name="poisson_2d.py-64"></a><span class="n">Gamma_left</span> <span class="o">=</span> <span class="n">LeftBoundary</span><span class="p">()</span>
<a name="poisson_2d.py-65"></a><span class="n">Gamma_left</span><span class="o">.</span><span class="n">mark</span><span class="p">(</span><span class="n">boundary_parts</span><span class="p">,</span> <span class="mi">2</span><span class="p">)</span>
<a name="poisson_2d.py-66"></a>
<a name="poisson_2d.py-67"></a><span class="c"># Mark right boundary as subdomain 3</span>
<a name="poisson_2d.py-68"></a><span class="k">class</span> <span class="nc">RightBoundary</span><span class="p">(</span><span class="n">SubDomain</span><span class="p">):</span>
<a name="poisson_2d.py-69"></a>    <span class="k">def</span> <span class="nf">inside</span><span class="p">(</span><span class="bp">self</span><span class="p">,</span> <span class="n">x</span><span class="p">,</span> <span class="n">on_boundary</span><span class="p">):</span>
<a name="poisson_2d.py-70"></a>        <span class="n">tol</span> <span class="o">=</span> <span class="mf">1E-16</span>   <span class="c"># tolerance for coordinate comparisons</span>
<a name="poisson_2d.py-71"></a>        <span class="k">return</span> <span class="n">on_boundary</span> <span class="ow">and</span> <span class="nb">abs</span><span class="p">(</span><span class="n">x</span><span class="p">[</span><span class="mi">0</span><span class="p">]</span> <span class="o">-</span> <span class="mi">1</span><span class="p">)</span> <span class="o">&lt;</span> <span class="n">tol</span>
<a name="poisson_2d.py-72"></a><span class="n">Gamma_right</span> <span class="o">=</span> <span class="n">RightBoundary</span><span class="p">()</span>
<a name="poisson_2d.py-73"></a><span class="n">Gamma_right</span><span class="o">.</span><span class="n">mark</span><span class="p">(</span><span class="n">boundary_parts</span><span class="p">,</span> <span class="mi">3</span><span class="p">)</span>
<a name="poisson_2d.py-74"></a>
<a name="poisson_2d.py-75"></a><span class="n">u_L</span> <span class="o">=</span> <span class="n">Expression</span><span class="p">(</span><span class="s">&#39;0&#39;</span><span class="p">)</span>
<a name="poisson_2d.py-76"></a><span class="n">u_R</span> <span class="o">=</span> <span class="n">Expression</span><span class="p">(</span><span class="s">&#39;0&#39;</span><span class="p">)</span>
<a name="poisson_2d.py-77"></a>
<a name="poisson_2d.py-78"></a><span class="n">sigma_bottom</span> <span class="o">=</span> <span class="n">Expression</span><span class="p">(</span><span class="s">&#39;-(pi/2.0)*sin(2*pi*x[0])&#39;</span><span class="p">)</span>
<a name="poisson_2d.py-79"></a><span class="n">sigma_top</span>    <span class="o">=</span> <span class="n">Expression</span><span class="p">(</span><span class="s">&#39;0&#39;</span><span class="p">)</span>
<a name="poisson_2d.py-80"></a>
<a name="poisson_2d.py-81"></a><span class="n">bcs</span> <span class="o">=</span> <span class="p">[</span><span class="n">DirichletBC</span><span class="p">(</span><span class="n">V</span><span class="p">,</span> <span class="n">u_L</span><span class="p">,</span> <span class="n">boundary_parts</span><span class="p">,</span> <span class="mi">2</span><span class="p">),</span>
<a name="poisson_2d.py-82"></a>       <span class="n">DirichletBC</span><span class="p">(</span><span class="n">V</span><span class="p">,</span> <span class="n">u_R</span><span class="p">,</span> <span class="n">boundary_parts</span><span class="p">,</span> <span class="mi">3</span><span class="p">)]</span>
<a name="poisson_2d.py-83"></a>
<a name="poisson_2d.py-84"></a><span class="n">u_e</span> <span class="o">=</span> <span class="n">Expression</span><span class="p">(</span><span class="s">&#39;sin(2*pi*x[0])*sin((pi/2.0)*x[1])&#39;</span><span class="p">)</span>
<a name="poisson_2d.py-85"></a>
<a name="poisson_2d.py-86"></a><span class="c"># Define variational problem</span>
<a name="poisson_2d.py-87"></a><span class="n">u</span> <span class="o">=</span> <span class="n">TrialFunction</span><span class="p">(</span><span class="n">V</span><span class="p">)</span>
<a name="poisson_2d.py-88"></a><span class="n">v</span> <span class="o">=</span> <span class="n">TestFunction</span><span class="p">(</span><span class="n">V</span><span class="p">)</span>
<a name="poisson_2d.py-89"></a><span class="n">f</span> <span class="o">=</span> <span class="n">Expression</span><span class="p">(</span><span class="s">&#39;(4.0*pi*pi+pi*pi/4.0)*(sin(2*pi*x[0])*sin((pi/2.0)*x[1]))&#39;</span><span class="p">)</span>
<a name="poisson_2d.py-90"></a><span class="n">a</span> <span class="o">=</span> <span class="n">inner</span><span class="p">(</span><span class="n">nabla_grad</span><span class="p">(</span><span class="n">u</span><span class="p">),</span> <span class="n">nabla_grad</span><span class="p">(</span><span class="n">v</span><span class="p">))</span><span class="o">*</span><span class="n">dx</span>
<a name="poisson_2d.py-91"></a><span class="n">L</span> <span class="o">=</span> <span class="n">f</span><span class="o">*</span><span class="n">v</span><span class="o">*</span><span class="n">dx</span> <span class="o">+</span> <span class="n">sigma_bottom</span><span class="o">*</span><span class="n">v</span><span class="o">*</span><span class="n">ds</span><span class="p">(</span><span class="mi">0</span><span class="p">)</span> <span class="o">+</span> <span class="n">sigma_top</span><span class="o">*</span><span class="n">v</span><span class="o">*</span><span class="n">ds</span><span class="p">(</span><span class="mi">1</span><span class="p">)</span>
<a name="poisson_2d.py-92"></a>
<a name="poisson_2d.py-93"></a><span class="n">M_equ</span> <span class="o">=</span> <span class="n">inner</span><span class="p">(</span><span class="n">u</span><span class="p">,</span> <span class="n">v</span><span class="p">)</span> <span class="o">*</span> <span class="n">dx</span>
<a name="poisson_2d.py-94"></a><span class="n">M</span> <span class="o">=</span> <span class="n">assemble</span><span class="p">(</span><span class="n">M_equ</span><span class="p">)</span>
<a name="poisson_2d.py-95"></a>
<a name="poisson_2d.py-96"></a><span class="c"># Compute solution</span>
<a name="poisson_2d.py-97"></a><span class="c">#A = assemble(a, exterior_facet_domains=boundary_parts)</span>
<a name="poisson_2d.py-98"></a><span class="c">#b = assemble(L, exterior_facet_domains=boundary_parts)</span>
<a name="poisson_2d.py-99"></a><span class="c">#for condition in bcs: condition.apply(A, b)</span>
<a name="poisson_2d.py-100"></a>
<a name="poisson_2d.py-101"></a><span class="n">A</span><span class="p">,</span> <span class="n">b</span> <span class="o">=</span> <span class="n">assemble_system</span><span class="p">(</span><span class="n">a</span><span class="p">,</span> <span class="n">L</span><span class="p">,</span> <span class="n">bcs</span><span class="p">,</span> <span class="n">exterior_facet_domains</span><span class="o">=</span><span class="n">boundary_parts</span><span class="p">)</span>
<a name="poisson_2d.py-102"></a><span class="k">if</span> <span class="n">mesh</span><span class="o">.</span><span class="n">num_cells</span><span class="p">()</span> <span class="o">&lt;</span> <span class="mi">64</span><span class="p">:</span> <span class="c"># print for small meshes only</span>
<a name="poisson_2d.py-103"></a>    <span class="k">print</span> <span class="n">A</span><span class="o">.</span><span class="n">array</span><span class="p">()</span>      <span class="c"># be careful with the identation!</span>
<a name="poisson_2d.py-104"></a>    <span class="k">print</span> <span class="n">M</span><span class="o">.</span><span class="n">array</span><span class="p">()</span>
<a name="poisson_2d.py-105"></a>    <span class="k">print</span> <span class="n">b</span><span class="o">.</span><span class="n">array</span><span class="p">()</span>
<a name="poisson_2d.py-106"></a>
<a name="poisson_2d.py-107"></a>
<a name="poisson_2d.py-108"></a>
<a name="poisson_2d.py-109"></a><span class="c"># Solve the linear system</span>
<a name="poisson_2d.py-110"></a><span class="n">u</span> <span class="o">=</span> <span class="n">Function</span><span class="p">(</span><span class="n">V</span><span class="p">)</span>
<a name="poisson_2d.py-111"></a><span class="c">#parameters[&quot;linear_algebra_backend&quot;] = &quot;PETSc&quot;</span>
<a name="poisson_2d.py-112"></a><span class="n">solver</span> <span class="o">=</span> <span class="n">KrylovSolver</span><span class="p">(</span><span class="s">&#39;cg&#39;</span><span class="p">,</span> <span class="s">&#39;ilu&#39;</span><span class="p">)</span>
<a name="poisson_2d.py-113"></a><span class="n">solver</span><span class="o">.</span><span class="n">parameters</span><span class="p">[</span><span class="s">&#39;absolute_tolerance&#39;</span><span class="p">]</span> <span class="o">=</span> <span class="mf">1E-14</span>
<a name="poisson_2d.py-114"></a><span class="n">solver</span><span class="o">.</span><span class="n">parameters</span><span class="p">[</span><span class="s">&#39;relative_tolerance&#39;</span><span class="p">]</span> <span class="o">=</span> <span class="mf">1E-14</span>
<a name="poisson_2d.py-115"></a><span class="n">solver</span><span class="o">.</span><span class="n">parameters</span><span class="p">[</span><span class="s">&#39;maximum_iterations&#39;</span><span class="p">]</span> <span class="o">=</span> <span class="mi">10000</span>
<a name="poisson_2d.py-116"></a><span class="n">set_log_level</span><span class="p">(</span><span class="n">PROGRESS</span><span class="p">)</span>
<a name="poisson_2d.py-117"></a>
<a name="poisson_2d.py-118"></a><span class="n">solver</span><span class="o">.</span><span class="n">solve</span><span class="p">(</span><span class="n">A</span><span class="p">,</span> <span class="n">u</span><span class="o">.</span><span class="n">vector</span><span class="p">(),</span> <span class="n">b</span><span class="p">)</span>
<a name="poisson_2d.py-119"></a>
<a name="poisson_2d.py-120"></a><span class="c"># Dump solution to file in VTK format</span>
<a name="poisson_2d.py-121"></a><span class="nb">file</span> <span class="o">=</span> <span class="n">File</span><span class="p">(</span><span class="s">&#39;poisson.pvd&#39;</span><span class="p">)</span>
<a name="poisson_2d.py-122"></a><span class="nb">file</span> <span class="o">&lt;&lt;</span> <span class="n">u</span>
<a name="poisson_2d.py-123"></a>
<a name="poisson_2d.py-124"></a><span class="n">plot</span><span class="p">(</span><span class="n">mesh</span><span class="p">)</span>
<a name="poisson_2d.py-125"></a><span class="n">plot</span><span class="p">(</span><span class="n">u</span><span class="p">)</span>
<a name="poisson_2d.py-126"></a><span class="n">interactive</span><span class="p">()</span>
</pre></div>
</td></tr></table>
    </div>
  


          </div>
          
        </div>
      </div>
    </div>
    <div data-modules="source/set-changeset" data-hash="909ad9e757141bc9ac07faec7c68a018b9a6a3c9"></div>




  
  
    <script id="branch-dialog-template" type="text/html">
  

<div class="tabbed-filter-widget branch-dialog">
  <div class="tabbed-filter">
    <input placeholder="Filter branches" class="filter-box" autosave="branch-dropdown-13334003" type="text">
    [[^ignoreTags]]
      <div class="aui-tabs horizontal-tabs aui-tabs-disabled filter-tabs">
        <ul class="tabs-menu">
          <li class="menu-item active-tab"><a href="#branches">Branches</a></li>
          <li class="menu-item"><a href="#tags">Tags</a></li>
        </ul>
      </div>
    [[/ignoreTags]]
  </div>
  
    <div class="tab-pane active-pane" id="branches" data-filter-placeholder="Filter branches">
      <ol class="filter-list">
        <li class="empty-msg">No matching branches</li>
        [[#branches]]
          
            [[#hasMultipleHeads]]
              [[#heads]]
                <li class="comprev filter-item">
                  <a class="pjax-trigger filter-item-link" href="/npetra/math292f15/src/[[changeset]]/codes/fenics/forward/poisson_2d.py?at=[[safeName]]"
                     title="[[name]]">
                    [[name]] ([[shortChangeset]])
                  </a>
                </li>
              [[/heads]]
            [[/hasMultipleHeads]]
            [[^hasMultipleHeads]]
              <li class="comprev filter-item">
                <a class="pjax-trigger filter-item-link" href="/npetra/math292f15/src/[[changeset]]/codes/fenics/forward/poisson_2d.py?at=[[safeName]]" title="[[name]]">
                  [[name]]
                </a>
              </li>
            [[/hasMultipleHeads]]
          
        [[/branches]]
      </ol>
    </div>
    <div class="tab-pane" id="tags" data-filter-placeholder="Filter tags">
      <ol class="filter-list">
        <li class="empty-msg">No matching tags</li>
        [[#tags]]
          <li class="comprev filter-item">
            <a class="pjax-trigger filter-item-link" href="/npetra/math292f15/src/[[changeset]]/codes/fenics/forward/poisson_2d.py?at=[[safeName]]" title="[[name]]">
              [[name]]
            </a>
          </li>
        [[/tags]]
      </ol>
    </div>
  
</div>

</script>
  



  </div>

        
        
        
          <script id="image-upload-template" type="text/html">
  

<form id="upload-image" method="POST"
    
      action="/xhr/npetra/math292f15/image-upload/"
    >
  <input type='hidden' name='csrfmiddlewaretoken' value='ONGo9v2kb2lAkOCXYybaFat0cK1pZ0xh' />

  <div class="drop-target">
    <p class="centered">Drag image here</p>
  </div>

  
  <div>
    <button class="aui-button click-target">Select an image</button>
    <input name="file" type="file" class="hidden file-target"
           accept="image/jpeg, image/gif, image/png" />
    <input type="submit" class="hidden">
  </div>
</form>


</script>
        
      </div>
    </div>
  </div>

    </div>
  </div>

  <footer id="footer" role="contentinfo" data-modules="components/footer">
    <section class="footer-body">
      
  <ul>
  <li>
    <a class="support-ga" target="_blank"
       data-support-gaq-page="Blog"
       href="http://blog.bitbucket.org">Blog</a>
  </li>
  <li>
    <a class="support-ga" target="_blank"
       data-support-gaq-page="Home"
       href="/support">Support</a>
  </li>
  <li>
    <a class="support-ga"
       data-support-gaq-page="PlansPricing"
       href="/plans">Plans &amp; pricing</a>
  </li>
  <li>
    <a class="support-ga" target="_blank"
       data-support-gaq-page="DocumentationHome"
       href="//confluence.atlassian.com/display/BITBUCKET">Documentation</a>
  </li>
  <li>
    <a class="support-ga" target="_blank"
       data-support-gaq-page="DocumentationAPI"
       href="//confluence.atlassian.com/x/IYBGDQ">API</a>
  </li>
  <li>
    <a class="support-ga" target="_blank"
       data-support-gaq-page="SiteStatus"
       href="http://status.bitbucket.org/">Site status</a>
  </li>
  <li>
    <a class="support-ga" id="meta-info"
       data-support-gaq-page="MetaInfo"
       href="#">Version info</a>
  </li>
  <li>
    <a class="support-ga" target="_blank"
       data-support-gaq-page="EndUserAgreement"
       href="//www.atlassian.com/end-user-agreement?utm_source=bitbucket&amp;utm_medium=link&amp;utm_campaign=footer">Terms of service</a>
  </li>
  <li>
    <a class="support-ga" target="_blank"
       data-support-gaq-page="PrivacyPolicy"
       href="//www.atlassian.com/company/privacy?utm_source=bitbucket&amp;utm_medium=link&amp;utm_campaign=footer">Privacy policy</a>
  </li>
</ul>
<div id="meta-info-content" style="display: none;">
  <ul>
    
      <li><a href="/account/user/uvilla/" class="view-language-link">English</a></li>
    
    <li>
      <a class="support-ga" target="_blank"
         data-support-gaq-page="GitDocumentation"
         href="http://git-scm.com/">Git 2.1.1</a>
    </li>
    <li>
      <a class="support-ga" target="_blank"
         data-support-gaq-page="HgDocumentation"
         href="http://mercurial.selenic.com/">Mercurial 2.9</a>
    </li>
    <li>
      <a class="support-ga" target="_blank"
         data-support-gaq-page="DjangoDocumentation"
         href="https://www.djangoproject.com/">Django 1.7.8</a>
    </li>
    <li>
      <a class="support-ga" target="_blank"
         data-support-gaq-page="PythonDocumentation"
         href="http://www.python.org/">Python 2.7.3</a>
    </li>
    <li>
      <a class="support-ga" target="_blank"
         data-support-gaq-page="DeployedVersion"
         href="#">b3d49699260a / b3d49699260a @ app-106</a>
    </li>
  </ul>
</div>
<ul class="atlassian-links">
  <li>
    <a id="atlassian-jira-link" target="_blank"
       title="Track everything – bugs, tasks, deadlines, code – and pull reports to stay informed."
       href="http://www.atlassian.com/software/jira?utm_source=bitbucket&amp;utm_medium=link&amp;utm_campaign=bitbucket_footer">JIRA</a>
  </li>
  <li>
    <a id="atlassian-confluence-link" target="_blank"
       title="Content Creation, Collaboration & Knowledge Sharing for Teams."
       href="http://www.atlassian.com/software/confluence/overview/team-collaboration-software?utm_source=bitbucket&amp;utm_medium=link&amp;utm_campaign=bitbucket_footer">Confluence</a>
  </li>
  <li>
    <a id="atlassian-bamboo-link" target="_blank"
       title="Continuous integration and deployment, release management."
       href="http://www.atlassian.com/software/bamboo/overview?utm_source=bitbucket&amp;utm_medium=link&amp;utm_campaign=bitbucket_footer">Bamboo</a>
  </li>
  <li>
    <a id="atlassian-sourcetree-link" target="_blank"
       title="A free Git and Mercurial desktop client for Mac or Windows."
       href="http://www.sourcetreeapp.com/?utm_source=bitbucket&amp;utm_medium=link&amp;utm_campaign=bitbucket_footer">SourceTree</a>
  </li>
  <li>
    <a id="atlassian-hipchat-link" target="_blank"
       title="Group chat and IM built for teams."
       href="http://www.hipchat.com/?utm_source=bitbucket&amp;utm_medium=link&amp;utm_campaign=bitbucket_footer">HipChat</a>
  </li>
</ul>
<div id="footer-logo">
  <a target="_blank"
     title="Bitbucket is developed by Atlassian in San Francisco and Austin."
     href="http://www.atlassian.com?utm_source=bitbucket&amp;utm_medium=logo&amp;utm_campaign=bitbucket_footer">Atlassian</a>
</div>

    </section>
  </footer>
</div>


  
  


<script src="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/jsi18n/en/djangojs.js"></script>
<script src="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/dist/webpack/vendor.js"></script>
<script src="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/dist/webpack/app.js"></script>


<script>
  (function () {
    var ga = document.createElement('script');
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    ga.setAttribute('async', 'true');
    var s = document.getElementsByTagName('script')[0];
    s.parentNode.insertBefore(ga, s);
  }());
</script>


  

<div data-modules="components/mentions/index">
  <script id="mention-result" type="text/html">
    
<div class="aui-avatar aui-avatar-small">
  <div class="aui-avatar-inner">
    <img src="[[avatar_url]]">
  </div>
</div>
[[#display_name]]
  <span class="display-name">[[&display_name]]</span> <small class="username">[[&username]]</small>
[[/display_name]]
[[^display_name]]
  <span class="username">[[&username]]</span>
[[/display_name]]
[[#is_teammate]][[^is_team]]
  <span class="aui-lozenge aui-lozenge-complete aui-lozenge-subtle">teammate</span>
[[/is_team]][[/is_teammate]]

  </script>
  <script id="mention-call-to-action" type="text/html">
    
[[^query]]
<li class="bb-typeahead-item">Begin typing to search for a user</li>
[[/query]]
[[#query]]
<li class="bb-typeahead-item">Continue typing to search for a user</li>
[[/query]]

  </script>
  <script id="mention-no-results" type="text/html">
    
[[^searching]]
<li class="bb-typeahead-item">Found no matching users for <em>[[query]]</em>.</li>
[[/searching]]
[[#searching]]
<li class="bb-typeahead-item bb-typeahead-searching">Searching for <em>[[query]]</em>.</li>
[[/searching]]

  </script>
</div>
<div data-modules="components/typeahead/emoji/index">
  <script id="emoji-result" type="text/html">
    
<div class="aui-avatar aui-avatar-small">
  <div class="aui-avatar-inner">
    <img src="[[src]]">
  </div>
</div>
<span class="name">[[&name]]</span>

  </script>
</div>

<div data-modules="components/repo-typeahead/index">
  <script id="repo-typeahead-result" type="text/html">
    <span class="aui-avatar aui-avatar-project aui-avatar-xsmall">
  <span class="aui-avatar-inner">
    <img src="[[avatar]]">
  </span>
</span>
<span class="owner">[[&owner]]</span>/<span class="slug">[[&slug]]</span>

  </script>
</div>
<script id="share-form-template" type="text/html">
    

<div class="error aui-message hidden">
  <span class="aui-icon icon-error"></span>
  <div class="message"></div>
</div>
<form class="aui">
  <table class="widget bb-list aui">
    <thead>
    <tr class="assistive">
      <th class="user">User</th>
      <th class="role">Role</th>
      <th class="actions">Actions</th>
    </tr>
    </thead>
    <tbody>
      <tr class="form">
        <td colspan="2">
          <input type="text" class="text bb-user-typeahead user-or-email"
                 placeholder="Username or email address"
                 autocomplete="off"
                 data-bb-typeahead-focus="false"
                 [[#disabled]]disabled[[/disabled]]>
        </td>
        <td class="actions">
          <button type="submit" class="aui-button aui-button-light" disabled>Add</button>
        </td>
      </tr>
    </tbody>
  </table>
</form>

  </script>
<script id="share-detail-template" type="text/html">
    

[[#username]]
<td class="user
    [[#hasCustomGroups]]custom-groups[[/hasCustomGroups]]"
    [[#error]]data-error="[[error]]"[[/error]]>
  <div title="[[displayName]]">
    <a href="/[[username]]/" class="user">
      <span class="aui-avatar aui-avatar-xsmall">
        <span class="aui-avatar-inner">
          <img src="[[avatar]]">
        </span>
      </span>
      <span>[[displayName]]</span>
    </a>
  </div>
</td>
[[/username]]
[[^username]]
<td class="email
    [[#hasCustomGroups]]custom-groups[[/hasCustomGroups]]"
    [[#error]]data-error="[[error]]"[[/error]]>
  <div title="[[email]]">
    <span class="aui-icon aui-icon-small aui-iconfont-email"></span>
    [[email]]
  </div>
</td>
[[/username]]
<td class="role
    [[#hasCustomGroups]]custom-groups[[/hasCustomGroups]]">
  <div>
    [[#group]]
      [[#hasCustomGroups]]
        <select class="group [[#noGroupChoices]]hidden[[/noGroupChoices]]">
          [[#groups]]
            <option value="[[slug]]"
                [[#isSelected]]selected[[/isSelected]]>
              [[name]]
            </option>
          [[/groups]]
        </select>
      [[/hasCustomGroups]]
      [[^hasCustomGroups]]
      <label>
        <input type="checkbox" class="admin"
            [[#isAdmin]]checked[[/isAdmin]]>
        Administrator
      </label>
      [[/hasCustomGroups]]
    [[/group]]
    [[^group]]
      <ul>
        <li class="permission aui-lozenge aui-lozenge-complete
            [[^read]]aui-lozenge-subtle[[/read]]"
            data-permission="read">
          read
        </li>
        <li class="permission aui-lozenge aui-lozenge-complete
            [[^write]]aui-lozenge-subtle[[/write]]"
            data-permission="write">
          write
        </li>
        <li class="permission aui-lozenge aui-lozenge-complete
            [[^admin]]aui-lozenge-subtle[[/admin]]"
            data-permission="admin">
          admin
        </li>
      </ul>
    [[/group]]
  </div>
</td>
<td class="actions
    [[#hasCustomGroups]]custom-groups[[/hasCustomGroups]]">
  <div>
    <a href="#" class="delete">
      <span class="aui-icon aui-icon-small aui-iconfont-remove">Delete</span>
    </a>
  </div>
</td>

  </script>
<script id="share-team-template" type="text/html">
    

<div class="clearfix">
  <span class="team-avatar-container">
    <span class="aui-avatar aui-avatar-medium">
      <span class="aui-avatar-inner">
        <img src="[[avatar]]">
      </span>
    </span>
  </span>
  <span class="team-name-container">
    [[display_name]]
  </span>
</div>
<p class="helptext">
  
    Existing users are granted access to this team immediately.
    New users will be sent an invitation.
  
</p>
<div class="share-form"></div>

  </script>
<script id="scope-list-template" type="text/html">
    <ul class="scope-list">
  [[#scopes]]
    <li class="scope-list--item">
      <span class="scope-list--icon aui-icon aui-icon-small [[icon]]"></span>
      [[description]]
    </li>
  [[/scopes]]
</ul>

  </script>


  

<script id="source-changeset" type="text/html">
  

<a href="/npetra/math292f15/src/[[raw_node]]/[[path]]?at=master"
    class="[[#selected]]highlight[[/selected]]"
    data-hash="[[node]]">
  [[#author.username]]
    <span class="aui-avatar aui-avatar-xsmall">
      <span class="aui-avatar-inner">
        <img src="[[author.avatar]]">
      </span>
    </span>
    <span class="author" title="[[raw_author]]">[[author.display_name]]</span>
  [[/author.username]]
  [[^author.username]]
    <span class="aui-avatar aui-avatar-xsmall">
      <span class="aui-avatar-inner">
        <img src="https://d3oaxc4q5k2d6q.cloudfront.net/m/b3d49699260a/img/default_avatar/16/user_blue.png">
      </span>
    </span>
    <span class="author unmapped" title="[[raw_author]]">[[author]]</span>
  [[/author.username]]
  <time datetime="[[utctimestamp]]" data-title="true">[[utctimestamp]]</time>
  <span class="message">[[message]]</span>
</a>

</script>
<script id="embed-template" type="text/html">
  

<form class="aui inline-dialog-embed-dialog">
  <label for="embed-code-[[dialogId]]">Embed this source in another page:</label>
  <input type="text" readonly="true" value="&lt;script src=&quot;[[url]]&quot;&gt;&lt;/script&gt;" id="embed-code-[[dialogId]]" class="embed-code">
</form>

</script>


  <script id="edit-form-template" type="text/html">
  


<form class="bb-content-container online-edit-form aui"
      data-repository="[[owner]]/[[slug]]"
      data-destination-repository="[[destinationOwner]]/[[destinationSlug]]"
      data-local-id="[[localID]]"
      data-is-writer="[[#isWriter]]true[[/isWriter]][[^isWriter]]false[[/isWriter]]"
      data-has-push-access="[[#hasPushAccess]]true[[/hasPushAccess]][[^hasPushAccess]]false[[/hasPushAccess]]"
      data-is-pull-request="[[#isPullRequest]]true[[/isPullRequest]][[^isPullRequest]]false[[/isPullRequest]]"
      data-hash="[[hash]]"
      data-branch="[[branch]]"
      data-path="[[path]]"
      data-is-create="[[isCreate]]"
      data-preview-url="/xhr/[[owner]]/[[slug]]/preview/[[hash]]/[[encodedPath]]"
      data-preview-error="We had trouble generating your preview."
      data-unsaved-changes-error="Your changes will be lost. Are you sure you want to leave?">
  <div class="bb-content-container-header">
    <div class="bb-content-container-header-primary">
      <h2 class="bb-content-container-heading">
        [[#isCreate]]
          
            Creating <span class="edit-path">[[path]]</span> on branch: <strong>[[branch]]</strong>
          
        [[/isCreate]]
        [[^isCreate]]
          
            Editing <span class="edit-path">[[path]]</span> on branch: <strong>[[branch]]</strong>
          
        [[/isCreate]]
      </h2>
    </div>
    <div class="bb-content-container-header-secondary">
      <div class="hunk-nav aui-buttons">
        <button class="prev-hunk-button aui-button aui-button-light"
              disabled="disabled" aria-disabled="true" title="previous change">&#x25B2;</button>
        <button class="next-hunk-button aui-button aui-button-light"
              disabled="disabled" aria-disabled="true" title="next change">&#x25BC;</button>
      </div>
    </div>
  </div>
  <div class="bb-content-container-body has-header has-footer file-editor">
    <textarea id="id_source"></textarea>
  </div>
  <div class="preview-pane"></div>
  <div class="bb-content-container-footer">
    <div class="bb-content-container-footer-primary">
      <div id="syntax-mode" class="bb-content-container-item field">
        <label for="id_syntax-mode">Syntax mode:</label>
        <select id="id_syntax-mode">
          [[#syntaxes]]
            <option value="[[#mime]][[mime]][[/mime]][[^mime]][[mode]][[/mime]]">[[name]]</option>
          [[/syntaxes]]
        </select>
      </div>
      <div id="indent-mode" class="bb-content-container-item field">
        <label for="id_indent-mode">Indent mode:</label>
        <select id="id_indent-mode">
          <option value="tabs">Tabs</option>
          <option value="spaces">Spaces</option>
        </select>
      </div>
      <div id="indent-size" class="bb-content-container-item field">
        <label for="id_indent-size">Indent size:</label>
        <select id="id_indent-size">
          <option value="2">2</option>
          <option value="4">4</option>
          <option value="8">8</option>
        </select>
      </div>
    </div>
    <div class="bb-content-container-footer-secondary">
      [[^isCreate]]
        <button class="preview-button aui-button aui-button-light"
                disabled="disabled" aria-disabled="true"
                data-preview-label="View diff"
                data-edit-label="Edit file">View diff</button>
      [[/isCreate]]
      <button class="save-button aui-button aui-button-primary"
              disabled="disabled" aria-disabled="true">Commit</button>
      [[^isCreate]]
        <a class="aui-button aui-button-link cancel-link" href="#">Cancel</a>
      [[/isCreate]]
    </div>
  </div>
</form>

</script>
  <script id="commit-form-template" type="text/html">
  

<form class="aui commit-form"
      data-title="Commit changes"
      [[#isDelete]]
        data-default-message="[[filename]] deleted online with Bitbucket"
      [[/isDelete]]
      [[#isCreate]]
        data-default-message="[[filename]] created online with Bitbucket"
      [[/isCreate]]
      [[^isDelete]]
        [[^isCreate]]
          data-default-message="[[filename]] edited online with Bitbucket"
        [[/isCreate]]
      [[/isDelete]]
      data-fork-error="We had trouble creating your fork."
      data-commit-error="We had trouble committing your changes."
      data-pull-request-error="We had trouble creating your pull request."
      data-update-error="We had trouble updating your pull request."
      data-branch-conflict-error="A branch with that name already exists."
      data-forking-message="Forking repository"
      data-committing-message="Committing changes"
      data-merging-message="Branching and merging changes"
      data-creating-pr-message="Creating pull request"
      data-updating-pr-message="Updating pull request"
      data-cta-label="Commit"
      data-cancel-label="Cancel">
  [[#isDelete]]
    <div class="aui-message info">
      <span class="aui-icon icon-info"></span>
      <span class="message">
        
          Committing this change will delete [[filename]] from your repository.
        
      </span>
    </div>
  [[/isDelete]]
  <div class="aui-message error hidden">
    <span class="aui-icon icon-error"></span>
    <span class="message"></span>
  </div>
  [[^isWriter]]
    <div class="aui-message info">
      <span class="aui-icon icon-info"></span>
      <p class="title">
        
          You don't have write access to this repository.
        
      </p>
      <span class="message">
        
          We'll create a fork for your changes and submit a
          pull request back to this repository.
        
      </span>
    </div>
  [[/isWriter]]
  [[#isRename]]
    <div class="field-group">
      <label for="id_path">New path</label>
      <input type="text" id="id_path" class="text" value="[[path]]"/>
    </div>
  [[/isRename]]
  <div class="field-group">
    <label for="id_message">Commit message</label>
    <textarea id="id_message" class="long-field textarea"></textarea>
  </div>
  [[^isPullRequest]]
    [[#isWriter]]
      <fieldset class="group">
        <div class="checkbox">
          [[#hasPushAccess]]
            <input id="id_create-pullrequest" class="checkbox" type="checkbox">
            <label for="id_create-pullrequest">Create a pull request for this change</label>
          [[/hasPushAccess]]
          [[^hasPushAccess]]
            <input id="id_create-pullrequest" class="checkbox" type="checkbox" checked="checked" aria-disabled="true" disabled="true">
            <label for="id_create-pullrequest" title="Branch restrictions do not allow you to update this branch.">Create a pull request for this change</label>
          [[/hasPushAccess]]

        </div>
      </fieldset>
      <div id="pr-fields">
        <div id="branch-name-group" class="field-group">
          <label for="id_branch-name">Branch name</label>
          <input type="text" id="id_branch-name" class="text long-field">
        </div>

        

<div class="field-group" id="id_reviewers_group">
  <label for="reviewers">Reviewers</label>

  
  <input id="reviewers" name="reviewers" type="hidden"
          value=""
          data-mention-url="/xhr/mentions/repositories/:dest_username/:dest_slug"
          data-reviewers="[]"
          data-suggested="[]"
          data-locked="[]">

  <div class="error"></div>
  <div class="suggested-reviewers"></div>

</div>


      </div>
    [[/isWriter]]
  [[/isPullRequest]]
  <button type="submit" id="id_submit">Commit</button>
</form>

</script>
  <script id="merge-message-template" type="text/html">
  Merged [[hash]] into [[branch]]

[[message]]

</script>
  <script id="commit-merge-error-template" type="text/html">
  



  We had trouble merging your changes. We stored them on the <strong>[[branch]]</strong> branch, so feel free to
  <a href="/[[owner]]/[[slug]]/full-commit/[[hash]]/[[path]]?at=[[encodedBranch]]">view them</a> or
  <a href="#" class="create-pull-request-link">create a pull request</a>.


</script>
  <script id="selected-reviewer-template" type="text/html">
  <div class="aui-avatar aui-avatar-xsmall">
  <div class="aui-avatar-inner">
    <img src="[[avatar_url]]">
  </div>
</div>
[[display_name]]

</script>
  <script id="suggested-reviewer-template" type="text/html">
  <button class="aui-button aui-button-link" type="button" tabindex="-1">[[display_name]]</button>

</script>
  <script id="suggested-reviewers-template" type="text/html">
  

<span class="suggested-reviewer-list-label">Recent:</span>
<ul class="suggested-reviewer-list unstyled-list"></ul>

</script>

  


  <div id="recently-mentioned-5236897"
       data-value="[]"></div>





  
  
  




<aui-inline-dialog
      id="super-touch-point-dialog"
      data-aui-alignment="bottom right"

      
      data-aui-alignment-static="true"
      data-modules="header/help-menu,connect/connect-views,connect/super-touch-point"
      responds-to="toggle"
      aria-hidden="true">

    <div id="ace-stp-section" class="no-touch-point">
      <div id="ace-stp-help-section">
        <h1 class="ace-stp-heading">Help</h1>

        <form id="ace-stp-search-form" class="aui" target="_blank" method="get"
            action="https://support.atlassian.com/customer/search">
          <span id="stp-search-icon" class="aui-icon aui-icon-large aui-iconfont-search"></span>
          <input id="ace-stp-search-form-input" name="q" class="text" type="text" placeholder="Ask a question">
        </form>

        <ul id="ace-stp-help-links">
          <li>
            <a class="support-ga" data-support-gaq-page="DocumentationHome"
                href="https://confluence.atlassian.com/x/bgozDQ" target="_blank">
              Online help
            </a>
          </li>
          <li>
            <a class="support-ga" data-support-gaq-page="GitTutorials"
                href="https://www.atlassian.com/git?utm_source=bitbucket&amp;utm_medium=link&amp;utm_campaign=help_dropdown&amp;utm_content=learn_git"
                target="_blank">
              Learn Git
            </a>
          </li>
          <li>
            <a id="keyboard-shortcuts-link"
               href="#">Keyboard shortcuts</a>
          </li>
          <li>
            <a href="/whats-new/" id="features-link">
              Latest features
            </a>
          </li>
          <li>
            <a class="support-ga" data-support-gaq-page="DocumentationTutorials"
                href="https://confluence.atlassian.com/x/Q4sFLQ" target="_blank">
              Bitbucket tutorials
            </a>
          </li>
          <li>
            <a class="support-ga" data-support-gaq-page="SiteStatus"
                href="http://status.bitbucket.org/" target="_blank">
              Site status
            </a>
          </li>
          <li>
            <a class="support-ga" data-support-gaq-page="Home" href="/support">
              Support
            </a>
          </li>

        </ul>
      </div>

      
    </div>
</aui-inline-dialog>
  


  <div class="omnibar" data-modules="components/omnibar/index">
    <form class="omnibar-form aui"></form>
  </div>
  <script id="omnibar-form-template" type="text/html">
    <div class="omnibar-input-container">
  <input class="omnibar-input" type="text">
</div>
<ul class="omnibar-result-group-list"></ul>

  </script>
  <script id="omnibar-blank-slate-template" type="text/html">
    

<div class="omnibar-blank-slate">No results found</div>

  </script>
  <script id="omnibar-result-group-list-item-template" type="text/html">
    <div class="omnibar-result-group-header clearfix">
  <h2 class="omnibar-result-group-label" title="[[label]]">[[label]]</h2>
  <span class="omnibar-result-group-context" title="[[context]]">[[context]]</span>
</div>
<ul class="omnibar-result-list unstyled-list"></ul>

  </script>
  <script id="omnibar-result-list-item-template" type="text/html">
    <span class="omnibar-result-label">[[&label]]</span>
[[#context]]
  <span class="omnibar-result-context">[[context]]</span>
[[/context]]

  </script>




<script type="text/javascript">window.NREUM||(NREUM={});NREUM.info={"beacon":"bam.nr-data.net","queueTime":0,"licenseKey":"a2cef8c3d3","agent":"js-agent.newrelic.com/nr-768.min.js","transactionName":"Z11RZxdWW0cEVkYLDV4XdUYLVEFdClsdAAtEWkZQDlJBGgRFQhFMQl1DXFcZQ10AQkFYBFlUVlEXWEJHAA==","userAttributes":"SxpaQDpWQEANUFwWC1NZR1YBFQ9AF0BXTkBZS2xSFV4XDgNUXhEHHBpGQABFaloEWFdAWBJNRVoJW1QWGA==","applicationID":"1841284","errorBeacon":"bam.nr-data.net","applicationTime":180}</script>
</body>
</html>
