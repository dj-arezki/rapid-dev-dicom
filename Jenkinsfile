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

      stage('make check') {
        steps {
          script {
            sh "make check"
          }

        }
      }

      stage('make dist') {
        steps {
          script {
            sh "make dist"
          }

        }
      }

      stage('Publish') {
        when {
          expression {
            return env.BRANCH_NAME ==~ /^development_.*/ || env.BRANCH_NAME ==~ /^release_.*/ || env.BRANCH_NAME == "main"
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
      options {
        disableConcurrentBuilds()
      }
    }