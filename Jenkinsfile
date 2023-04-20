/*******************************************************************************
 * IBM Confidential
 *
 * (C) Copyright IBM Corp. 2019, 2020, 2021 All Rights Reserved.
 *
 * US Government Users Restricted Rights - Use, duplication or
 * disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
 *
 *******************************************************************************/
String whi_64gb_node = 'whi-test-01'
String whi_scan_agent = 'whis-bld02'
String PRJ_NAME = 'cal-core'
String DOCKER_REPOSITORY = 'wh-imaging-dev-docker-local.artifactory.swg-devops.com'
BUILD_UTIL_IMAGE = "whi-algo-build-util:latest"

// AppScan on Cloud vars
String APPSCAN_APP_ID = '6e257cb6-52da-4eca-a4f8-c81900886793'

pipeline {

    agent{
        label whi_scan_agent
    }

    options{
        disableConcurrentBuilds()
    }
    stages{
        stage('Prep') {
            steps {
                withCredentials([
                    usernamePassword(credentialsId:'WHIDevOps-jfrog-cred',
                                     passwordVariable: 'WHI_REPO_PASSWORD', usernameVariable:'WHI_REPO_USER')]){
                    script {
                        deleteDir()
                        checkout scm
                        sh "git clean -dfx"
                    }
                }
            }
        }

        stage('make'){
            steps {
                withCredentials([usernamePassword(credentialsId: 'WHIDevOps-jfrog-cred', passwordVariable: 'WHI_REPO_PASSWORD', usernameVariable: 'WHI_REPO_USER')]) {
                    script {
                        sh "make"
                    }
                }
            }
        }

        stage('make check'){
            steps {
                script {
                    sh "make check"
                }
            }
        }

        stage('make dist'){
            steps {
                script {
                    sh "make dist"
                }
            }
        }

        stage('Publish'){
            when { expression { return env.BRANCH_NAME ==~ /^development_.*/ || env.BRANCH_NAME ==~ /^release_.*/ || env.BRANCH_NAME == "main" } }
            steps {
                withCredentials([
                    usernamePassword(credentialsId: 'WHIDevOps-jfrog-cred', passwordVariable: 'WHI_REPO_PASSWORD', usernameVariable: 'WHI_REPO_USER')]) {
                    script {
                        sh "make publish"
                    }
                }
            }
        }

    }
}
